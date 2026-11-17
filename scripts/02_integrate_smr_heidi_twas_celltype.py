#!/usr/bin/env python3
"""Real Multi-Evidence Causal Triangulation for Scheme 2 (SMR/HEIDI + TWAS + MAGMA + snRNA-seq).

Integrates strictly real summary outputs staged under ../inputs/:
  - TWAS_union_87_exported_candidates.tsv  (TWAS Z/P + MAGMA + BrainMeta SMR/HEIDI across 87 genes)
  - TWAS_12_shared_candidates_cross_method.tsv (12 cross-method shared candidates)
  - SMR_DCC_SLC12A5_ANKK1_NEK4.tsv         (BrainMeta strict SMR/HEIDI metrics)
  - coloc_ready_loci.tsv                   (DNAH10, SYNE1, RAB2A)
  - gene_level_cross_trait_all.tsv.gz      (PGC MDD & Anxiety MAGMA)
  - pseudobulk_de_results.csv              (DLPFC snRNA-seq MDD differential expression across 11 cell types)

Separates causal candidates (Bonferroni SMR + HEIDI P >= 0.05) from linkage-confounded
associations (HEIDI P < 0.05) and maps their DLPFC cell-type specificity.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent


def classify_causal_evidence(df: pd.DataFrame) -> pd.DataFrame:
    """Apply strict SMR + HEIDI + cross-trait tiering rules."""
    smr_p = pd.to_numeric(df["MDD_SMR_P"], errors="coerce")
    heidi_p = pd.to_numeric(df["MDD_HEIDI"], errors="coerce")
    anx_smr_p = pd.to_numeric(df["ANX_SMR_P"], errors="coerce")
    anx_heidi_p = pd.to_numeric(df["ANX_HEIDI"], errors="coerce")

    tiers = []
    for i in range(len(df)):
        sp_m = smr_p.iloc[i]
        hp_m = heidi_p.iloc[i]
        sp_a = anx_smr_p.iloc[i]
        hp_a = anx_heidi_p.iloc[i]
        shared_smr = str(df["Shared_SMR_primary"].iloc[i]).lower() == "true"

        if shared_smr and pd.notna(hp_m) and hp_m >= 0.05 and pd.notna(hp_a) and hp_a >= 0.05:
            tiers.append("Tier1_Shared_Causal_SMR_HEIDI_Pass")
        elif pd.notna(sp_m) and sp_m < 5e-5 and pd.notna(hp_m) and hp_m >= 0.05:
            tiers.append("Tier2_MDD_Causal_SMR_HEIDI_Pass")
        elif pd.notna(sp_m) and sp_m < 5e-5 and pd.notna(hp_m) and hp_m < 0.05:
            tiers.append("Tier3_Linkage_Confounded_HEIDI_Fail")
        else:
            tiers.append("Tier4_TWAS_Support_Only")
    df = df.copy()
    df["Causal_Evidence_Tier"] = tiers
    return df


def summarize_snrna_celltypes(pb: pd.DataFrame, target_genes: set[str]) -> pd.DataFrame:
    """For each target gene, extract top MDD-perturbed cell type and mean expression level."""
    sub = pb[pb["gene"].isin(target_genes)].copy()
    sub["abs_stat"] = sub["stat"].abs()
    top_cell = (
        sub.sort_values("pvalue")
        .drop_duplicates("gene")
        .rename(columns={
            "cell_type": "top_snRNA_DE_celltype",
            "log2FoldChange": "top_snRNA_log2FC",
            "stat": "top_snRNA_stat",
            "pvalue": "top_snRNA_p",
            "padj": "top_snRNA_padj",
            "baseMean": "top_snRNA_baseMean",
        })[["gene", "top_snRNA_DE_celltype", "top_snRNA_log2FC", "top_snRNA_stat", "top_snRNA_p", "top_snRNA_padj", "top_snRNA_baseMean"]]
    )
    highest_expr = (
        sub.sort_values("baseMean", ascending=False)
        .drop_duplicates("gene")[["gene", "cell_type", "baseMean"]]
        .rename(columns={"cell_type": "highest_expr_celltype", "baseMean": "highest_expr_baseMean"})
    )
    return top_cell.merge(highest_expr, on="gene", how="left"), sub


def plot_causal_triangulation(tri_df: pd.DataFrame, sub_pb: pd.DataFrame, fig_dir: Path) -> None:
    """Plot SMR vs HEIDI causality filter and DLPFC cell-type expression heatmap for top candidates."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # Panel A: SMR -log10(P) vs HEIDI P showing causal quadrant vs linkage quadrant
    df_p = tri_df.dropna(subset=["MDD_SMR_P", "MDD_HEIDI"]).copy()
    df_p["neglog10_smr"] = -np.log10(df_p["MDD_SMR_P"].astype(float).clip(lower=1e-30))
    df_p["heidi_val"] = df_p["MDD_HEIDI"].astype(float)

    tier_palette = {
        "Tier1_Shared_Causal_SMR_HEIDI_Pass": ("#187C92", "D", 60, "Tier 1: Shared SMR + HEIDI Pass"),
        "Tier2_MDD_Causal_SMR_HEIDI_Pass": ("#29688B", "o", 38, "Tier 2: MDD SMR + HEIDI Pass"),
        "Tier3_Linkage_Confounded_HEIDI_Fail": ("#C66B45", "X", 48, "Tier 3: HEIDI Failed (Linkage)"),
        "Tier4_TWAS_Support_Only": ("#8D9AA2", "o", 22, "Tier 4: TWAS Only"),
    }

    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.2), dpi=300)

    ax0 = axes[0]
    for tier, (col, mk, sz, lbl) in tier_palette.items():
        sub = df_p[df_p["Causal_Evidence_Tier"] == tier]
        if len(sub):
            ax0.scatter(sub["heidi_val"], sub["neglog10_smr"], c=col, marker=mk, s=sz, alpha=0.88, edgecolor="#233746", linewidth=0.5, label=lbl)

    ax0.axvline(0.05, color="#B83B2A", linestyle="--", linewidth=1.1, label="HEIDI P = 0.05 cutoff")
    ax0.axhline(-np.log10(5e-5), color="#667782", linestyle=":", linewidth=1.0)
    # Annotate key validated genes (Tier 1 + notable Tier 3 negative controls like ANKK1/NEK4)
    offsets = {
        "DCC": (-26, 6),
        "GPX1": (7, -3),
        "NEK4": (5, 4),
        "ANKK1": (5, -2),
    }
    ann_genes = {"SLC12A5", "DCC", "AMT", "GPX1", "PPP6C", "ANKK1", "NEK4", "SPCS1", "RNF123"}
    for _, row in df_p[df_p["Gene"].isin(ann_genes)].iterrows():
        ox, oy = offsets.get(row["Gene"], (5, 3))
        ax0.annotate(row["Gene"], (row["heidi_val"], row["neglog10_smr"]), xytext=(ox, oy), textcoords="offset points", fontsize=7.8, fontweight="bold", color="#233746")

    ax0.set_xlabel(r"HEIDI Test $P$-value ($\geq 0.05$ supports shared pleiotropy/causality)")
    ax0.set_ylabel(r"BrainMeta MDD $-\log_{10}(P_{\mathrm{SMR}})$")
    ax0.set_title("A. BrainMeta SMR Causality vs. HEIDI Linkage Filter", fontweight="bold", loc="left", color="#233746")
    ax0.legend(frameon=False, fontsize=7.5, loc="upper right")

    # Panel B: Heatmap of normalized snRNA pseudobulk Wald z-stat across 11 cortical cell types for prioritized causal & coloc genes
    focus_genes = ["SLC12A5", "DCC", "PPP6C", "GPX1", "SPCS1", "RNF123", "SCAMP1", "GNL3", "SYNE1", "DNAH10", "RAB2A", "ANKK1"]
    hm = sub_pb[sub_pb["gene"].isin(focus_genes)].pivot_table(index="gene", columns="cell_type", values="stat", aggfunc="first")
    hm = hm.reindex([g for g in focus_genes if g in hm.index])
    ax1 = axes[1]
    im = ax1.imshow(hm.fillna(0.0).to_numpy(), cmap="RdBu_r", vmin=-3.0, vmax=3.0, aspect="auto")
    ax1.set_xticks(np.arange(len(hm.columns)))
    ax1.set_xticklabels(hm.columns, rotation=45, ha="right", fontsize=8)
    ax1.set_yticks(np.arange(len(hm.index)))
    ax1.set_yticklabels(hm.index, fontsize=8.5, fontweight="bold")
    for r_i in range(len(hm.index)):
        for c_i in range(len(hm.columns)):
            val = hm.iloc[r_i, c_i]
            txt = f"{val:.1f}" if pd.notna(val) else "NA"
            ax1.text(c_i, r_i, txt, ha="center", va="center", fontsize=6.8, color="#111111" if abs(val if pd.notna(val) else 0) < 2.0 else "white")
    cbar = fig.colorbar(im, ax=ax1, fraction=0.045, pad=0.04)
    cbar.set_label("snRNA MDD Pseudobulk Wald stat")
    ax1.set_title("B. DLPFC Cortical Cell-Type MDD Perturbation of Causal Candidates", fontweight="bold", loc="left", color="#233746")

    fig.tight_layout(pad=1.5)
    fig.savefig(fig_dir / "Fig2_Scheme2_SMR_HEIDI_and_CellType_Map.png", dpi=300)
    plt.close(fig)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="SMR/HEIDI + TWAS + MAGMA + DLPFC snRNA causal triangulation.")
    parser.add_argument("--in-dir", type=Path, default=DEFAULT_BASE_DIR / "inputs")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_BASE_DIR / "outputs")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    args = parse_args(argv)
    in_dir: Path = args.in_dir
    out_dir: Path = args.out_dir
    tab_dir = out_dir / "tables"
    fig_dir = out_dir / "figures"
    bm_dir = out_dir / "benchmarks"
    for d in (tab_dir, fig_dir, bm_dir):
        d.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    twas_87 = pd.read_csv(in_dir / "TWAS_union_87_exported_candidates.tsv", sep="\t")
    coloc_top3 = pd.read_csv(in_dir / "coloc_ready_loci.tsv", sep="\t")
    smr_detail = pd.read_csv(in_dir / "SMR_DCC_SLC12A5_ANKK1_NEK4.tsv", sep="\t")
    magma = pd.read_csv(in_dir / "gene_level_cross_trait_all.tsv.gz", sep="\t")
    pb = pd.read_csv(in_dir / "pseudobulk_de_results.csv")

    tri = classify_causal_evidence(twas_87)
    extra_smr_rows = []
    for _, r in smr_detail.iterrows():
        if r["symbol"] not in set(tri["Gene"]):
            extra_smr_rows.append({
                "Gene": r["symbol"],
                "Chr": 18 if r["symbol"] == "DCC" else np.nan,
                "MAGMA_MDD_P": r["FUMA_P_MDD"],
                "MAGMA_ANX_P": r["FUMA_P_ANX"],
                "MAGMA_both": str(r["FUMA_both_Bonf"]).lower() == "true",
                "MDD_SMR_P": r["MDD_p_SMR"],
                "ANX_SMR_P": r["ANX_p_SMR"],
                "MDD_HEIDI": r["MDD_p_HEIDI"],
                "ANX_HEIDI": r["ANX_p_HEIDI"],
                "Shared_SMR_primary": str(r["shared_primary"]).lower() == "true",
                "Causal_Evidence_Tier": (
                    "Tier1_Shared_Causal_SMR_HEIDI_Pass"
                    if (str(r["shared_primary"]).lower() == "true" and r["MDD_p_HEIDI"] >= 0.05 and r["ANX_p_HEIDI"] >= 0.05)
                    else "Tier2_MDD_Causal_SMR_HEIDI_Pass"
                ),
            })
    extra_smr_df = pd.DataFrame(extra_smr_rows)

    known_genes = set(tri["Gene"]) | (set(extra_smr_df["Gene"]) if len(extra_smr_df) else set())
    extra_coloc = coloc_top3[~coloc_top3["gene"].isin(known_genes)][["gene", "chr", "lead_snp_id", "lead_snp_p"]].rename(
        columns={"gene": "Gene", "chr": "Chr"}
    )
    extra_coloc["Causal_Evidence_Tier"] = "Coloc_Top3_FineMapped_Locus"
    combined = pd.concat([tri, extra_smr_df, extra_coloc], ignore_index=True)

    mag_clean = magma.dropna(subset=["SYMBOL"]).drop_duplicates("SYMBOL")[
        ["SYMBOL", "ZSTAT_MDD", "P_MDD", "ZSTAT_ANX", "P_ANX", "Evidence_tier"]
    ]
    combined = combined.merge(mag_clean, left_on="Gene", right_on="SYMBOL", how="left")

    target_set = set(combined["Gene"].dropna())
    snrna_sum, sub_pb = summarize_snrna_celltypes(pb, target_set)
    combined = combined.merge(snrna_sum, left_on="Gene", right_on="gene", how="left")
    combined["snRNA_provenance_note"] = "May2026_Leiden_r0.8_unadjusted_sex_pending_raw_GSE213982_41subtype_rerun"

    combined = combined.sort_values(["Causal_Evidence_Tier", "MDD_SMR_P"], na_position="last").reset_index(drop=True)
    combined.to_csv(tab_dir / "causal_triangulation_master_table.tsv", sep="\t", index=False)

    tier_counts = combined["Causal_Evidence_Tier"].value_counts().reset_index()
    tier_counts.columns = ["Causal_Evidence_Tier", "n_genes"]
    tier_counts.to_csv(tab_dir / "causal_tier_counts.tsv", sep="\t", index=False)

    if len(sub_pb):
        plot_causal_triangulation(combined, sub_pb, fig_dir)

    elapsed = time.perf_counter() - t0
    with open(bm_dir / "triangulation_benchmark.json", "w") as f:
        json.dump({"elapsed_seconds": round(elapsed, 3), "n_genes_triangulated": int(len(combined))}, f, indent=2)

    print("Scheme 2 Step 2 (SMR/HEIDI + TWAS + MAGMA + snRNA Triangulation) finished.")
    print(tier_counts.to_string(index=False))
    return combined, tier_counts


if __name__ == "__main__":
    main()
