#!/usr/bin/env python3
"""Local Downstream Integration Pipeline for Scheme 1 (snRNA-seq + Spatial Calibration + PGC MAGMA).

Integrates staged local datasets under ../inputs/:
  1. pseudobulk_de_results.csv            (259,817 DLPFC snRNA-seq pseudobulk DE tests across 11 clusters)
  2. spatial_deconvolution_and_layers.csv (3,639 Visium spots from May 2026 single-donor 151673 benchmark)
  3. gene_level_cross_trait_all.tsv.gz    (17,451 PGC MDD & Anxiety MAGMA gene-level statistics)
  4. OXPHOS_GO0006119_148genes.tsv        (Unified GO:0006119 canonical symbol set from v3.22.0)
  5. MitoCarta_core_frozen.tsv            (100 nuclear MitoCarta3.0 OXPHOS structural subunits, optional)
  6. SMR_DCC_SLC12A5_ANKK1_NEK4.tsv & TWAS_12_shared_candidates_cross_method.tsv

Important scientific provenance notes enforced explicitly in outputs:
  - May 2026 pseudobulk_de_results.csv used Leiden r=0.8 cluster mapping (~diagnosis without sex covariate);
    cluster 'Pericyte' contains excitatory neuronal markers (e.g., TBR1). Reprocessing raw GSE213982
    with embedded 8 major / 41 author subtype barcode tokens is scripted in the Colab Pro notebook.
  - May 2026 spatial_deconvolution_and_layers.csv originates from ONE healthy neurotypical donor
    (spatialLIBD 151673) with synthetic random spot-level diagnosis labels and equal y-bin cuts;
    spatial coordinates/counts provide an authentic Moran's I spatial-autocorrelation calibration
    benchmark, NOT patient-control MDD spatial pathology.
  - Cell-type OXPHOS MWU/Wilcoxon P-values assume gene independence (analogous to fixed-correlation
    CAMERA rho=0.01) and are reported alongside single-gene FDR<0.05 counts (= 0 across all types).
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import scipy.sparse as sp
import scipy.stats as stats
import statsmodels.api as sm
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests


DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent


def current_rss_mb() -> float:
    """Return current process resident set size in MiB."""
    return psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)


def compute_morans_i(coords: np.ndarray, values: np.ndarray, k: int = 6) -> tuple[float, float, float]:
    """Compute k-NN row-standardized Moran's I and analytical z-score using exact sparse weights."""
    n = len(values)
    values = np.asarray(values, dtype=float)
    z = values - values.mean()
    denom = float(np.sum(z ** 2))
    if denom == 0.0 or n < 5:
        return 0.0, 0.0, 1.0

    k_eff = min(k, n - 1)
    nbrs = NearestNeighbors(n_neighbors=k_eff + 1).fit(coords)
    _, indices = nbrs.kneighbors(coords)
    nbr_idx = indices[:, 1 : k_eff + 1]

    row_ind = np.repeat(np.arange(n, dtype=np.int32), k_eff)
    col_ind = nbr_idx.ravel().astype(np.int32)
    data = np.full(n * k_eff, 1.0 / k_eff, dtype=float)
    W = sp.csr_matrix((data, (row_ind, col_ind)), shape=(n, n))

    num = float(z @ (W @ z))
    s0 = float(W.sum())
    s1 = 0.5 * float((W + W.T).power(2).sum())
    row_sums = np.asarray(W.sum(axis=1)).ravel()
    col_sums = np.asarray(W.sum(axis=0)).ravel()
    s2 = float(np.sum((row_sums + col_sums) ** 2))

    moran_i = (n / s0) * (num / denom)
    ei = -1.0 / (n - 1)
    b2 = (n * float(np.sum(z ** 4))) / (denom ** 2)
    num_v = n * ((n ** 2 - 3 * n + 3) * s1 - n * s2 + 3 * (s0 ** 2)) - b2 * (
        (n ** 2 - n) * s1 - 2 * n * s2 + 6 * (s0 ** 2)
    )
    den_v = (n - 1) * (n - 2) * (n - 3) * (s0 ** 2)
    var_i = max((num_v / den_v) - (ei ** 2), 1e-12)
    z_score = (moran_i - ei) / np.sqrt(var_i)
    p_val = float(2.0 * (1.0 - stats.norm.cdf(abs(z_score))))
    return float(moran_i), float(z_score), p_val


def run_celltype_magma_integration(
    pb: pd.DataFrame,
    magma: pd.DataFrame,
    smr: pd.DataFrame,
    twas: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test whether cell-type MDD pseudobulk DE genes carry elevated MAGMA genetic liability."""
    mag_clean = (
        magma.dropna(subset=["SYMBOL", "ZSTAT_MDD", "Z_Stouffer_exploratory"])
        .sort_values("P_MDD")
        .drop_duplicates(subset=["SYMBOL"])
    )
    smr_tier1 = set(smr[smr["shared_primary"].astype(str).str.lower() == "true"]["symbol"])
    twas_12 = set(twas["Gene"].dropna())

    merged = pb.merge(mag_clean, left_on="gene", right_on="SYMBOL", how="inner")
    merged["abs_stat"] = merged["stat"].abs()
    merged["is_deg_fdr05"] = merged["padj"].fillna(1.0) < 0.05
    merged["is_deg_nom05"] = merged["pvalue"].fillna(1.0) < 0.05
    merged["is_pleio_bonf"] = merged["Evidence_tier"].eq("Bonferroni in both")

    enrich_rows = []
    overlap_frames = []

    for ctype, sub in merged.groupby("cell_type"):
        n_genes = len(sub)
        n_fdr05 = int(sub["is_deg_fdr05"].sum())
        n_nom05 = int(sub["is_deg_nom05"].sum())

        rho_mdd, p_rho_mdd = stats.spearmanr(sub["abs_stat"], sub["ZSTAT_MDD"])
        rho_stouf, p_rho_stouf = stats.spearmanr(sub["abs_stat"], sub["Z_Stouffer_exploratory"])

        q90 = sub["abs_stat"].quantile(0.90)
        top_grp = sub[sub["abs_stat"] >= q90]
        bg_grp = sub[sub["abs_stat"] < q90]
        _, p_top_mdd = stats.mannwhitneyu(
            top_grp["ZSTAT_MDD"], bg_grp["ZSTAT_MDD"], alternative="greater"
        )
        _, p_top_stouf = stats.mannwhitneyu(
            top_grp["Z_Stouffer_exploratory"], bg_grp["Z_Stouffer_exploratory"], alternative="greater"
        )

        a = int((sub["is_deg_nom05"] & sub["is_pleio_bonf"]).sum())
        b = int((sub["is_deg_nom05"] & ~sub["is_pleio_bonf"]).sum())
        c = int((~sub["is_deg_nom05"] & sub["is_pleio_bonf"]).sum())
        d = int((~sub["is_deg_nom05"] & ~sub["is_pleio_bonf"]).sum())
        or_pleio, p_fisher_pleio = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        p_fdr05_mdd = np.nan
        mean_z_fdr05 = np.nan
        if n_fdr05 >= 10:
            deg_sub = sub[sub["is_deg_fdr05"]]
            non_sub = sub[~sub["is_deg_fdr05"]]
            if len(non_sub) > 0:
                _, p_fdr05_mdd = stats.mannwhitneyu(
                    deg_sub["ZSTAT_MDD"], non_sub["ZSTAT_MDD"], alternative="greater"
                )
            mean_z_fdr05 = float(deg_sub["ZSTAT_MDD"].mean())

        cluster_note = (
            "Leiden_r0.8_cluster14_contains_TBR1_verify_on_raw_GSE213982"
            if str(ctype) == "Pericyte"
            else "May2026_Leiden_r0.8_unadjusted_sex_pending_raw_GSE213982_41subtype_rerun"
        )

        enrich_rows.append({
            "cell_type": ctype,
            "n_overlap_genes": n_genes,
            "n_deg_fdr05": n_fdr05,
            "n_deg_nom05": n_nom05,
            "spearman_abs_stat_vs_MAGMA_MDD": float(rho_mdd),
            "spearman_p_MDD": float(p_rho_mdd),
            "spearman_abs_stat_vs_Stouffer_MDD_ANX": float(rho_stouf),
            "spearman_p_Stouffer": float(p_rho_stouf),
            "top10pct_mean_MAGMA_Z_MDD": float(top_grp["ZSTAT_MDD"].mean()),
            "bg90pct_mean_MAGMA_Z_MDD": float(bg_grp["ZSTAT_MDD"].mean()),
            "top10pct_MWU_p_MDD": float(p_top_mdd),
            "top10pct_MWU_p_Stouffer": float(p_top_stouf),
            "deg_fdr05_mean_MAGMA_Z_MDD": mean_z_fdr05,
            "deg_fdr05_MWU_p_MDD": p_fdr05_mdd,
            "nomDE_x_BonfPleio_hits": a,
            "nomDE_x_BonfPleio_OR": float(or_pleio),
            "nomDE_x_BonfPleio_fisher_p": float(p_fisher_pleio),
            "upstream_annotation_caveat": cluster_note,
        })

        sub_ov = sub[
            (sub["is_deg_fdr05"] & (sub["MDD_nominal"].astype(str).str.lower() == "true"))
            | sub["gene"].isin(twas_12 | smr_tier1)
        ].copy()
        if len(sub_ov):
            sub_ov["in_TWAS12"] = sub_ov["gene"].isin(twas_12)
            sub_ov["in_SMR_Tier1"] = sub_ov["gene"].isin(smr_tier1)
            overlap_frames.append(
                sub_ov[[
                    "cell_type", "gene", "log2FoldChange", "stat", "pvalue", "padj",
                    "ZSTAT_MDD", "P_MDD", "Z_Stouffer_exploratory", "Evidence_tier",
                    "in_TWAS12", "in_SMR_Tier1",
                ]]
            )

    res_df = pd.DataFrame(enrich_rows).sort_values("top10pct_MWU_p_MDD").reset_index(drop=True)
    res_df["spearman_fdr_q_MDD"] = multipletests(res_df["spearman_p_MDD"], method="fdr_bh")[1]
    res_df["top10pct_fdr_q_MDD"] = multipletests(res_df["top10pct_MWU_p_MDD"], method="fdr_bh")[1]
    res_df["top10pct_fdr_q_Stouffer"] = multipletests(res_df["top10pct_MWU_p_Stouffer"], method="fdr_bh")[1]
    res_df["pleio_fisher_fdr_q"] = multipletests(res_df["nomDE_x_BonfPleio_fisher_p"], method="fdr_bh")[1]

    if overlap_frames:
        ov_df = pd.concat(overlap_frames, ignore_index=True).sort_values(["padj", "P_MDD"], na_position="last")
    else:
        ov_df = pd.DataFrame(columns=[
            "cell_type", "gene", "log2FoldChange", "stat", "pvalue", "padj",
            "ZSTAT_MDD", "P_MDD", "Z_Stouffer_exploratory", "Evidence_tier",
            "in_TWAS12", "in_SMR_Tier1",
        ])
    return res_df, ov_df


def run_celltype_oxphos_follow_up(
    pb: pd.DataFrame,
    oxphos_df: pd.DataFrame,
    mitocarta_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Audit GO:0006119 and MitoCarta core OXPHOS shifts across DLPFC snRNA clusters with independence caveats."""
    ox_genes = set(oxphos_df["Gene"].astype(str).str.strip())
    mito_genes = (
        set(mitocarta_df["MitoCartaSymbol"].astype(str).str.strip().str.replace('"', "", regex=False))
        if mitocarta_df is not None and "MitoCartaSymbol" in mitocarta_df.columns
        else set()
    )
    rows = []
    for ctype, sub in pb.groupby("cell_type"):
        sub_c = sub.dropna(subset=["stat", "log2FoldChange", "pvalue"]).copy()
        in_ox = sub_c["gene"].isin(ox_genes)
        ox_sub = sub_c[in_ox]
        bg_sub = sub_c[~in_ox]
        if len(ox_sub) < 5 or len(bg_sub) < 5:
            continue

        med_stat = float(ox_sub["stat"].median())
        mean_stat = float(ox_sub["stat"].mean())
        med_lfc = float(ox_sub["log2FoldChange"].median())
        pct_down = float((ox_sub["stat"] < 0).mean() * 100.0)
        _, p_wilcox = stats.wilcoxon(ox_sub["stat"])
        _, p_mwu_signed = stats.mannwhitneyu(ox_sub["stat"], bg_sub["stat"], alternative="two-sided")
        n_ox_fdr05 = int((ox_sub["padj"].fillna(1.0) < 0.05).sum())
        n_ox_nom05 = int((ox_sub["pvalue"].fillna(1.0) < 0.05).sum())

        mito_sub = sub_c[sub_c["gene"].isin(mito_genes)] if mito_genes else pd.DataFrame()
        mito_med_stat = float(mito_sub["stat"].median()) if len(mito_sub) >= 5 else np.nan
        mito_fdr05 = int((mito_sub["padj"].fillna(1.0) < 0.05).sum()) if len(mito_sub) else 0

        rows.append({
            "cell_type": ctype,
            "n_oxphos_covered": int(len(ox_sub)),
            "median_log2FC": med_lfc,
            "median_wald_stat": med_stat,
            "mean_wald_stat": mean_stat,
            "pct_downregulated": pct_down,
            "oxphos_fdr05_degs": n_ox_fdr05,
            "oxphos_nom05_degs": n_ox_nom05,
            "wilcoxon_signed_p_naive": float(p_wilcox),
            "mwu_vs_bg_p_naive": float(p_mwu_signed),
            "mitocarta100_covered": int(len(mito_sub)),
            "mitocarta100_median_wald": mito_med_stat,
            "mitocarta100_fdr05_degs": mito_fdr05,
            "stat_caveat": "Gene-level MWU/Wilcoxon omits intergene correlation across co-regulated OXPHOS genes; requires donor-level CAMERA/scDRS calibration",
        })
    df = pd.DataFrame(rows).sort_values("mwu_vs_bg_p_naive").reset_index(drop=True)
    if len(df):
        df["mwu_fdr_q_naive"] = multipletests(df["mwu_vs_bg_p_naive"], method="fdr_bh")[1]
        df["wilcoxon_fdr_q_naive"] = multipletests(df["wilcoxon_signed_p_naive"], method="fdr_bh")[1]
    return df


def run_spatial_layer_analysis(sp: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute spatial autocorrelation & layer calibration statistics on staged Visium table.

    Provenance guardrail: staged Visium file comes from May 2026 single-donor spatialLIBD 151673
    with synthetic per-spot diagnosis labels and equal-bin y-coordinate layer cuts.
    """
    provenance_tag = "calibration_benchmark_healthy_151673_synthetic_diagnosis_and_ybin_layers"
    layer_order = [l for l in ["L1", "L2", "L3", "L4", "L5", "L6", "WM"] if l in set(sp["cortical_layer"])]
    ct = pd.crosstab(sp["cortical_layer"], sp["diagnosis"]).reindex(layer_order, fill_value=0)
    chi2, p_chi2, _, _ = stats.chi2_contingency(ct) if ct.shape[0] >= 2 and ct.shape[1] == 2 else (0.0, 1.0, 0, None)

    tot_mdd = int(ct["MDD"].sum()) if "MDD" in ct.columns else 0
    tot_ctrl = int(ct["Control"].sum()) if "Control" in ct.columns else 0

    layer_rows = []
    for lyr in layer_order:
        m = int(ct.loc[lyr, "MDD"]) if "MDD" in ct.columns else 0
        c = int(ct.loc[lyr, "Control"]) if "Control" in ct.columns else 0
        or_val, p_fish = stats.fisher_exact([[m, max(0, tot_mdd - m)], [c, max(0, tot_ctrl - c)]])

        sub_l = sp[sp["cortical_layer"] == lyr].copy()
        mdd_g = sub_l[sub_l["diagnosis"] == "MDD"]["log1p_n_genes_by_counts"]
        ctrl_g = sub_l[sub_l["diagnosis"] == "Control"]["log1p_n_genes_by_counts"]
        if len(mdd_g) > 0 and len(ctrl_g) > 0:
            _, p_genes = stats.mannwhitneyu(mdd_g, ctrl_g, alternative="two-sided")
        else:
            p_genes = 1.0

        sub_l["is_mdd"] = (sub_l["diagnosis"] == "MDD").astype(float)
        if "_scvi_batch" in sub_l.columns and sub_l["_scvi_batch"].nunique() > 1:
            X = sm.add_constant(pd.get_dummies(sub_l[["is_mdd", "_scvi_batch"]], columns=["_scvi_batch"], drop_first=True, dtype=float))
        else:
            X = sm.add_constant(sub_l[["is_mdd"]], has_constant="add")
        ols = sm.OLS(sub_l["log1p_total_counts"].astype(float), X).fit(cov_type="HC3")

        layer_rows.append({
            "cortical_layer": lyr,
            "n_MDD_spots_synthetic": m,
            "n_Control_spots_synthetic": c,
            "pct_MDD_of_diagnosis": round(100.0 * m / max(1, tot_mdd), 2),
            "pct_Control_of_diagnosis": round(100.0 * c / max(1, tot_ctrl), 2),
            "composition_OR": float(or_val),
            "composition_fisher_p": float(p_fish),
            "MDD_mean_log1p_genes": float(mdd_g.mean()) if len(mdd_g) else np.nan,
            "Ctrl_mean_log1p_genes": float(ctrl_g.mean()) if len(ctrl_g) else np.nan,
            "log1p_genes_MWU_p": float(p_genes),
            "batch_adj_log1p_counts_beta": float(ols.params.get("is_mdd", np.nan)),
            "batch_adj_log1p_counts_p": float(ols.pvalues.get("is_mdd", 1.0)),
            "global_chi2_p": float(p_chi2),
            "provenance_status": provenance_tag,
        })
    layer_df = pd.DataFrame(layer_rows)
    if len(layer_df):
        layer_df["composition_fdr_q"] = multipletests(layer_df["composition_fisher_p"], method="fdr_bh")[1]
        layer_df["counts_ols_fdr_q"] = multipletests(layer_df["batch_adj_log1p_counts_p"], method="fdr_bh")[1]

    moran_rows = []
    # Overall healthy section Moran's I across all 3,639 real spots + stratified calibration check
    coords_all = sp[["array_row", "array_col"]].to_numpy(dtype=float)
    for metric in ["log1p_n_genes_by_counts", "log1p_total_counts"]:
        vals_all = sp[metric].to_numpy(dtype=float)
        mi, z_mi, p_mi = compute_morans_i(coords_all, vals_all, k=6)
        moran_rows.append({
            "stratum": "All_spots_healthy_151673",
            "n_spots": len(sp),
            "metric": metric,
            "k_neighbors": 6,
            "Moran_I": mi,
            "Moran_z": z_mi,
            "Moran_p": p_mi,
            "provenance_status": "authentic_spatialLIBD_151673_geometry",
        })
    for diag, sub_d in sp.groupby("diagnosis"):
        coords = sub_d[["array_row", "array_col"]].to_numpy(dtype=float)
        for metric in ["log1p_n_genes_by_counts", "log1p_total_counts"]:
            vals = sub_d[metric].to_numpy(dtype=float)
            mi, z_mi, p_mi = compute_morans_i(coords, vals, k=6)
            moran_rows.append({
                "stratum": f"SyntheticLabel_{diag}",
                "n_spots": len(sub_d),
                "metric": metric,
                "k_neighbors": 6,
                "Moran_I": mi,
                "Moran_z": z_mi,
                "Moran_p": p_mi,
                "provenance_status": provenance_tag,
            })
    return layer_df, pd.DataFrame(moran_rows)


def plot_scheme1_summary(
    cell_df: pd.DataFrame,
    ox_df: pd.DataFrame,
    layer_df: pd.DataFrame,
    fig_dir: Path,
) -> None:
    """Render publication-grade multi-panel summary plots for Scheme 1."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(1, 3, figsize=(16.8, 5.6), dpi=300)

    c_plot = cell_df.sort_values("top10pct_MWU_p_MDD", ascending=True)
    y = np.arange(len(c_plot))
    neglogp = -np.log10(c_plot["top10pct_MWU_p_MDD"].clip(lower=1e-300))
    colors = ["#C66B45" if q < 0.05 else "#187C92" if p < 0.05 else "#7D8B94"
              for p, q in zip(c_plot["top10pct_MWU_p_MDD"], c_plot["top10pct_fdr_q_MDD"])]
    axes[0].barh(y, neglogp, color=colors, edgecolor="#233746", linewidth=0.5)
    axes[0].axvline(-np.log10(0.05), color="#B83B2A", linestyle="--", linewidth=1.1, label="P = 0.05")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(c_plot["cell_type"])
    axes[0].invert_yaxis()
    axes[0].set_xlabel(r"$-\log_{10}(P)$ (Top 10% |DE| vs BG MAGMA $Z_{\mathrm{MDD}}$)")
    axes[0].set_title("A. DLPFC Cell-Type Genetic Enrichment", fontweight="bold", fontsize=10.5, loc="left", color="#233746")
    axes[0].legend(frameon=False, fontsize=8)

    x_pos = np.arange(len(c_plot))
    bars = axes[1].bar(x_pos, c_plot["n_deg_fdr05"], color="#29688B", alpha=0.85, edgecolor="#233746", linewidth=0.5)
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels([f"{c}*" if c == "Pericyte" else c for c in c_plot["cell_type"]], rotation=45, ha="right", fontsize=8)
    axes[1].set_ylabel("FDR < 0.05 MDD Pseudobulk DEGs (n)")
    axes[1].set_title("B. May 2026 Leiden DE (*Audit Caveat)", fontweight="bold", fontsize=10.5, loc="left", color="#233746")
    for b, n_v in zip(bars, c_plot["n_deg_fdr05"]):
        if n_v > 0:
            axes[1].text(b.get_x() + b.get_width() / 2, n_v + 4, str(int(n_v)), ha="center", va="bottom", fontsize=7.5)

    # Panel C: Honest Spatial Calibration Benchmark (Single Healthy Donor 151673; Synthetic Labels)
    ly_y = np.arange(len(layer_df))
    ors = layer_df["composition_OR"].to_numpy()
    l_cols = ["#7D8B94" if q >= 0.05 else "#C66B45" for q in layer_df["composition_fdr_q"]]
    axes[2].barh(ly_y, np.log2(ors), color=l_cols, edgecolor="#233746", linewidth=0.5)
    axes[2].axvline(0, color="#233746", linestyle="--", linewidth=1.0)
    axes[2].set_yticks(ly_y)
    axes[2].set_yticklabels(layer_df["cortical_layer"])
    axes[2].invert_yaxis()
    axes[2].set_xlabel(r"$\log_2(\mathrm{OR})$ (Synthetic Labels on Healthy 151673)")
    axes[2].set_title("C. Spatial Null (Healthy 151673; All q > 0.10)", fontweight="bold", fontsize=10.5, loc="left", color="#233746")

    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.16, wspace=0.28)
    fig.savefig(fig_dir / "Fig1_Scheme1_CellType_MAGMA_Spatial_Summary.png", dpi=300)
    plt.close(fig)

    if len(ox_df):
        ox_p = ox_df.sort_values("median_wald_stat")
        fig, ax = plt.subplots(figsize=(8.6, 5.1), dpi=300)
        yp = np.arange(len(ox_p))
        ax.barh(yp, ox_p["median_wald_stat"], color="#667782", edgecolor="#233746", linewidth=0.6)
        ax.axvline(0, color="#233746", linestyle="--", linewidth=1.0)
        ax.set_yticks(yp)
        ax.set_yticklabels([f"{ct} (n={n}, FDR<0.05: {d})" for ct, n, d in zip(ox_p["cell_type"], ox_p["n_oxphos_covered"], ox_p["oxphos_fdr05_degs"])])
        ax.set_xlabel("Median Pseudobulk Wald Statistic (GO:0006119 OXPHOS, MDD vs. Control)")
        ax.set_title("Canonical 148-Gene OXPHOS Across DLPFC Clusters\n(0 Single-Gene FDR < 0.05; Intergene Correlation Caveat)", fontweight="bold", fontsize=10.5, loc="center", color="#233746")
        fig.subplots_adjust(left=0.28, right=0.97, top=0.86, bottom=0.13)
        fig.savefig(fig_dir / "Fig2_Scheme1_Cortical_OXPHOS_by_CellType.png", dpi=300)
        plt.close(fig)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scheme 1 local snRNA + MAGMA + spatial calibration integration.")
    parser.add_argument("--in-dir", type=Path, default=DEFAULT_BASE_DIR / "inputs")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_BASE_DIR / "outputs")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> dict:
    args = parse_args(argv)
    in_dir: Path = args.in_dir
    out_dir: Path = args.out_dir
    tab_dir = out_dir / "tables"
    fig_dir = out_dir / "figures"
    bm_dir = out_dir / "benchmarks"
    for d in (tab_dir, fig_dir, bm_dir):
        d.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    rss_start = current_rss_mb()

    pb = pd.read_csv(in_dir / "pseudobulk_de_results.csv")
    sp_df = pd.read_csv(in_dir / "spatial_deconvolution_and_layers.csv", index_col=0)
    magma = pd.read_csv(in_dir / "gene_level_cross_trait_all.tsv.gz", sep="\t")
    oxphos = pd.read_csv(in_dir / "OXPHOS_GO0006119_148genes.tsv", sep="\t")
    smr = pd.read_csv(in_dir / "SMR_DCC_SLC12A5_ANKK1_NEK4.tsv", sep="\t")
    twas = pd.read_csv(in_dir / "TWAS_12_shared_candidates_cross_method.tsv", sep="\t")
    mito_path = in_dir / "MitoCarta_core_frozen.tsv"
    mitocarta = pd.read_csv(mito_path, sep="\t") if mito_path.exists() else None

    cell_df, ov_df = run_celltype_magma_integration(pb, magma, smr, twas)
    ox_df = run_celltype_oxphos_follow_up(pb, oxphos, mitocarta)
    layer_df, moran_df = run_spatial_layer_analysis(sp_df)
    plot_scheme1_summary(cell_df, ox_df, layer_df, fig_dir)

    cell_df.to_csv(tab_dir / "celltype_magma_enrichment_11types.tsv", sep="\t", index=False)
    ov_df.to_csv(tab_dir / "celltype_deg_genetics_overlap_genes.tsv", sep="\t", index=False)
    ox_df.to_csv(tab_dir / "celltype_oxphos_enrichment_11types.tsv", sep="\t", index=False)
    layer_df.to_csv(tab_dir / "spatial_layer_diagnosis_stats.tsv", sep="\t", index=False)
    moran_df.to_csv(tab_dir / "spatial_autocorrelation_moran.tsv", sep="\t", index=False)

    elapsed = time.perf_counter() - t0
    rss_end = current_rss_mb()
    meta = {
        "script": "01_local_celltype_magma_spatial_integration.py",
        "elapsed_seconds": round(elapsed, 3),
        "rss_start_mb": round(rss_start, 2),
        "rss_end_mb": round(rss_end, 2),
        "n_pseudobulk_rows": int(len(pb)),
        "n_spatial_spots": int(len(sp_df)),
        "n_magma_genes": int(len(magma)),
        "mitocarta_core_loaded": bool(mitocarta is not None),
    }
    with open(bm_dir / "downstream_execution_metrics.json", "w") as f:
        json.dump(meta, f, indent=2)

    print("Scheme 1 local integration completed cleanly.")
    print(json.dumps(meta, indent=2))
    return meta


if __name__ == "__main__":
    main()
