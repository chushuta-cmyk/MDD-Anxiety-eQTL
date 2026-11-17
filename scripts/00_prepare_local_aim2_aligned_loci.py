#!/usr/bin/env python3
"""
Aim 2 Local Preprocessing & Dual-Trait Large-N Locus Preparation
=================================================================
Runs completely locally on 16 GB RAM (< 1 GB peak RSS, < 10 s):
  1. Streams paired N=1.15M MDD (`MDD_EUR_GWAShug.tsv.gz`) and allele-aligned
     N=1.10M Anxiety (`ANX_EUR_GWAShug_aligned3.tsv.gz`) summary statistics.
  2. Extracts 7 expanded Aim 2 multi-signal/causal GRCh37 windows:
       - SLC12A5_chr20 (20:44.15-45.20 Mb)
       - DCC_chr18 (18:49.60-51.50 Mb)
       - PPP6C_chr9 (9:127.40-128.45 Mb)
       - SLC25A12_chr2 (2:172.39-173.12 Mb)
       - DNAH10_CDK2AP1_chr12 (12:123.00-124.60 Mb)
       - chr3p21_31_GPX1_AMT_RNF123 (3:48.70-50.10 Mb)
       - chr3p21_1_GNL3_NEK4_SPCS1 (3:52.10-53.30 Mb)
  3. Verifies A1/A2 alignment across MDD & ANX, computes per-trait Wakefield ABF
     95% credible sets under large-N GWAS power, and maps all 83 PP.H4 > 0.8
     COLOC leads from `coloc_leads.txt` / `coloc_aligned_out.txt` to CHR:BP.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import psutil

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR
REPO_ROOT = SCRIPT_DIR.parent
S2_ARCH = REPO_ROOT / "02_Scheme2_GWAS_Causal_Coloc_SMR" / "inputs" / "large_n_gwas_archive"

LOCUS_WINDOWS: dict[str, tuple[int, int, int]] = {
    "SLC12A5_chr20": (20, 44_150_000, 45_200_000),
    "DCC_chr18": (18, 49_600_000, 51_500_000),
    "PPP6C_chr9": (9, 127_400_000, 128_450_000),
    "SLC25A12_chr2": (2, 172_390_000, 173_120_000),
    "DNAH10_CDK2AP1_chr12": (12, 123_000_000, 124_600_000),
    "chr3p21_31_GPX1_AMT_RNF123": (3, 48_700_000, 50_100_000),
    "chr3p21_1_GNL3_NEK4_SPCS1": (3, 52_100_000, 53_300_000),
}


def wakefield_abf(beta: np.ndarray, se: np.ndarray, W: float = 0.04, snp_prior: float = 1e-4) -> tuple[np.ndarray, np.ndarray]:
    v = np.maximum(se ** 2, 1e-12)
    z = beta / np.sqrt(v)
    r = W / (W + v)
    log_abf = 0.5 * (np.log(1.0 - r) + r * (z ** 2))
    max_la = float(np.max(log_abf))
    unnorm = snp_prior * np.exp(log_abf - max_la)
    denom = float(np.exp(-max_la) + np.sum(unnorm))
    pip = unnorm / np.maximum(denom, 1e-300)
    # Posterior conditional on single causal variant (normalized across locus SNPs)
    pip_cond = np.exp(log_abf - max_la)
    pip_cond = pip_cond / float(np.sum(pip_cond))
    return pip_cond, np.asarray(pip, dtype=float)


def cs95_size(pip: np.ndarray) -> int:
    order = np.argsort(pip)[::-1]
    csum = np.cumsum(pip[order])
    hits = np.where(csum >= 0.95)[0]
    return int(hits[0] + 1) if len(hits) else int(len(pip))


def stream_gwas(
    gwas_path: Path,
    lead_set: set[str],
    windows: dict[str, tuple[int, int, int]],
    chunksize: int = 400_000,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, int]:
    """Stream gzipped/plain GWAShug table extracting both window slices and lead SNPs."""
    win_hits: dict[str, list[pd.DataFrame]] = {k: [] for k in windows}
    lead_hits: list[pd.DataFrame] = []
    n_total = 0

    usecols = ["CHR", "BP", "SNP", "A1", "A2", "BETA", "SE", "Z", "P", "N"]
    dtypes = {"SNP": "string", "A1": "string", "A2": "string"}
    for chunk in pd.read_csv(
        gwas_path,
        sep=r"\s+",
        usecols=lambda c: c in usecols,
        dtype=dtypes,
        chunksize=chunksize,
        low_memory=True,
    ):
        n_total += len(chunk)
        chunk["CHR"] = pd.to_numeric(chunk["CHR"], errors="coerce")
        chunk["BP"] = pd.to_numeric(chunk["BP"], errors="coerce")
        chunk = chunk.dropna(subset=["CHR", "BP", "SNP", "BETA", "SE", "P"])
        chunk["CHR"] = chunk["CHR"].astype(np.int32)
        chunk["BP"] = chunk["BP"].astype(np.int64)

        for tag, (c_target, start_bp, end_bp) in windows.items():
            sub = chunk[(chunk["CHR"] == c_target) & (chunk["BP"] >= start_bp) & (chunk["BP"] <= end_bp)]
            if len(sub):
                win_hits[tag].append(sub)

        if lead_set:
            lsub = chunk[chunk["SNP"].isin(lead_set)]
            if len(lsub):
                lead_hits.append(lsub)

    win_dfs = {
        tag: (pd.concat(parts, ignore_index=True).drop_duplicates("SNP") if parts else pd.DataFrame())
        for tag, parts in win_hits.items()
    }
    lead_df = pd.concat(lead_hits, ignore_index=True).drop_duplicates("SNP") if lead_hits else pd.DataFrame()
    return win_dfs, lead_df, n_total


def main(argv: list[str] | None = None) -> dict:
    p = argparse.ArgumentParser(description="Prepare paired large-N MDD+ANX dual-trait Aim 2 locus tables locally.")
    p.add_argument("--gwas-mdd", type=Path, default=S2_ARCH / "MDD_EUR_GWAShug.tsv.gz")
    p.add_argument("--gwas-anx", type=Path, default=S2_ARCH / "ANX_EUR_GWAShug_aligned3.tsv.gz")
    p.add_argument("--coloc-leads", type=Path, default=S2_ARCH / "coloc_leads.txt")
    p.add_argument("--coloc-out", type=Path, default=S2_ARCH / "coloc_aligned_out.txt")
    p.add_argument("--out-dir", type=Path, default=BASE_DIR / "outputs")
    p.add_argument("--chunksize", type=int, default=400_000)
    args = p.parse_args(argv)

    t0 = time.perf_counter()
    proc = psutil.Process()
    rss_start = proc.memory_info().rss / 1024**2

    loci_dir = args.out_dir / "loci_dualtrait"
    tab_dir = args.out_dir / "tables"
    bm_dir = args.out_dir / "benchmarks"
    for d in (loci_dir, tab_dir, bm_dir):
        d.mkdir(parents=True, exist_ok=True)

    lead_list: list[str] = []
    if args.coloc_leads.exists():
        lead_list = [ln.strip() for ln in args.coloc_leads.read_text().splitlines() if ln.strip()]
    lead_set = set(lead_list)

    mdd_wins, mdd_leads, n_mdd_snps = stream_gwas(args.gwas_mdd, lead_set, LOCUS_WINDOWS, chunksize=args.chunksize)
    rss_mid = proc.memory_info().rss / 1024**2
    anx_wins, anx_leads, n_anx_snps = stream_gwas(args.gwas_anx, lead_set, LOCUS_WINDOWS, chunksize=args.chunksize)
    rss_end = proc.memory_info().rss / 1024**2

    locus_summary_rows = []
    for tag, (chr_id, start_bp, end_bp) in LOCUS_WINDOWS.items():
        dm = mdd_wins.get(tag, pd.DataFrame())
        da = anx_wins.get(tag, pd.DataFrame())
        if dm.empty or da.empty:
            continue

        merged = dm.merge(
            da,
            on=["CHR", "BP", "SNP"],
            suffixes=("_MDD", "_ANX"),
            how="inner",
        )
        merged["A1_MDD_u"] = merged["A1_MDD"].str.upper()
        merged["A2_MDD_u"] = merged["A2_MDD"].str.upper()
        merged["A1_ANX_u"] = merged["A1_ANX"].str.upper()
        merged["A2_ANX_u"] = merged["A2_ANX"].str.upper()

        same = (merged["A1_MDD_u"] == merged["A1_ANX_u"]) & (merged["A2_MDD_u"] == merged["A2_ANX_u"])
        flip = (merged["A1_MDD_u"] == merged["A2_ANX_u"]) & (merged["A2_MDD_u"] == merged["A1_ANX_u"])
        n_direct = int(same.sum())
        n_flipped = int(flip.sum())

        merged = merged[same | flip].copy()
        flip_mask = (merged["A1_MDD_u"] == merged["A2_ANX_u"]).to_numpy()
        if flip_mask.any():
            merged.loc[flip_mask, "BETA_ANX"] = -merged.loc[flip_mask, "BETA_ANX"]
            merged.loc[flip_mask, "Z_ANX"] = -merged.loc[flip_mask, "Z_ANX"]

        merged = merged[(merged["SE_MDD"] > 0) & (merged["SE_ANX"] > 0)].sort_values("BP").reset_index(drop=True)
        if len(merged) < 10:
            continue

        # Recompute canonical signed Z = BETA / SE so stale Z columns in pre-aligned archives cannot invert sign
        stale_z_mask = (merged["Z_ANX"].to_numpy(float) * merged["BETA_ANX"].to_numpy(float)) < 0
        n_stale_z = int(stale_z_mask.sum())
        merged["Z_MDD"] = merged["BETA_MDD"] / merged["SE_MDD"]
        merged["Z_ANX"] = merged["BETA_ANX"] / merged["SE_ANX"]

        pip_mdd, _ = wakefield_abf(merged["BETA_MDD"].to_numpy(float), merged["SE_MDD"].to_numpy(float))
        pip_anx, _ = wakefield_abf(merged["BETA_ANX"].to_numpy(float), merged["SE_ANX"].to_numpy(float))
        merged["PIP_MDD_W004"] = pip_mdd
        merged["PIP_ANX_W004"] = pip_anx
        merged = merged.drop(columns=["A1_MDD_u", "A2_MDD_u", "A1_ANX_u", "A2_ANX_u"])
        merged.to_csv(loci_dir / f"{tag}_dualtrait_grch37.tsv", sep="\t", index=False)

        top_m_idx = int(np.argmin(merged["P_MDD"].to_numpy(float)))
        top_a_idx = int(np.argmin(merged["P_ANX"].to_numpy(float)))
        z_corr = float(np.corrcoef(merged["Z_MDD"].to_numpy(float), merged["Z_ANX"].to_numpy(float))[0, 1])

        locus_summary_rows.append({
            "locus_tag": tag,
            "chr": chr_id,
            "win_start": start_bp,
            "win_end": end_bp,
            "n_shared_snps": len(merged),
            "n_direct_allele_match": n_direct,
            "n_flipped_allele": n_flipped,
            "n_anx_stale_z_repaired": n_stale_z,
            "mdd_lead_snp": merged.loc[top_m_idx, "SNP"],
            "mdd_min_p": float(merged.loc[top_m_idx, "P_MDD"]),
            "mdd_lead_pip": float(merged.loc[top_m_idx, "PIP_MDD_W004"]),
            "mdd_cs95_size": cs95_size(pip_mdd),
            "anx_lead_snp": merged.loc[top_a_idx, "SNP"],
            "anx_min_p": float(merged.loc[top_a_idx, "P_ANX"]),
            "anx_lead_pip": float(merged.loc[top_a_idx, "PIP_ANX_W004"]),
            "anx_cs95_size": cs95_size(pip_anx),
            "local_z_pearson_r": round(z_corr, 4),
            "same_lead_snp": bool(merged.loc[top_m_idx, "SNP"] == merged.loc[top_a_idx, "SNP"]),
        })

    sum_df = pd.DataFrame(locus_summary_rows)
    sum_df.to_csv(tab_dir / "aim2_largeN_7loci_dualtrait_finemap.tsv", sep="\t", index=False)

    # Annotate 83 COLOC leads
    coloc_map = {}
    if args.coloc_out.exists():
        cout = pd.read_csv(args.coloc_out, sep="\t")
        for _, r in cout.iterrows():
            for snp_tok in str(r.get("CI_snplist", "")).split(","):
                s = snp_tok.strip()
                if s in lead_set and s not in coloc_map:
                    coloc_map[s] = (str(r["LOC"]), float(r["PP.H4.abf"]), int(r["nsnps"]))

    anx_leads_ren = anx_leads.rename(columns={
        "A1": "A1_ANX", "A2": "A2_ANX", "BETA": "BETA_ANX", "SE": "SE_ANX", "Z": "Z_ANX", "P": "P_ANX", "N": "N_ANX"
    })
    mdd_leads_ren = mdd_leads.rename(columns={
        "A1": "A1_MDD", "A2": "A2_MDD", "BETA": "BETA_MDD", "SE": "SE_MDD", "Z": "Z_MDD", "P": "P_MDD", "N": "N_MDD"
    })
    leads_df = pd.DataFrame({"SNP": lead_list})
    if not mdd_leads_ren.empty:
        leads_df = leads_df.merge(mdd_leads_ren, on="SNP", how="left")
    if not anx_leads_ren.empty:
        anx_cols = ["SNP", "A1_ANX", "A2_ANX", "BETA_ANX", "SE_ANX", "Z_ANX", "P_ANX", "N_ANX"]
        leads_df = leads_df.merge(anx_leads_ren[[c for c in anx_cols if c in anx_leads_ren.columns]], on="SNP", how="left")

    leads_df["coloc_LOC"] = leads_df["SNP"].map(lambda s: coloc_map.get(s, ("", np.nan, 0))[0])
    leads_df["coloc_PP_H4"] = leads_df["SNP"].map(lambda s: coloc_map.get(s, ("", np.nan, 0))[1])
    leads_df["coloc_nsnps"] = leads_df["SNP"].map(lambda s: coloc_map.get(s, ("", np.nan, 0))[2])
    if {"BETA_MDD", "BETA_ANX", "SE_MDD", "SE_ANX", "A1_MDD", "A1_ANX"}.issubset(leads_df.columns):
        same_l = leads_df["A1_MDD"].astype(str).str.upper() == leads_df["A1_ANX"].astype(str).str.upper()
        aligned_b = np.where(same_l, leads_df["BETA_ANX"], -leads_df["BETA_ANX"])
        leads_df["BETA_ANX_aligned"] = aligned_b
        leads_df["Z_ANX_aligned"] = aligned_b / leads_df["SE_ANX"]
        leads_df["same_effect_direction"] = (leads_df["BETA_MDD"] * aligned_b) > 0
    leads_df.to_csv(tab_dir / "aim2_83_coloc_leads_annotated.tsv", sep="\t", index=False)

    bench = {
        "script": "00_prepare_local_aim2_aligned_loci.py",
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "mdd_snps_streamed": int(n_mdd_snps),
        "anx_snps_streamed": int(n_anx_snps),
        "n_loci_extracted": int(len(sum_df)),
        "n_coloc_leads_matched_coords": int(leads_df["CHR"].notna().sum()) if "CHR" in leads_df.columns else 0,
        "rss_start_mb": round(rss_start, 2),
        "rss_mid_mb": round(rss_mid, 2),
        "rss_end_mb": round(rss_end, 2),
    }
    with open(bm_dir / "aim2_local_prep_benchmark.json", "w") as f:
        json.dump(bench, f, indent=2)

    print(json.dumps(bench, indent=2))
    return bench


if __name__ == "__main__":
    main()
