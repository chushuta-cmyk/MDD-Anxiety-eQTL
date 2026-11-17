#!/usr/bin/env python3
"""Real GWAS Locus Extraction & Single-Variant Wakefield ABF Fine-Mapping (Scheme 2, Step 1).

Streams the full 463.56 MB GRCh37 PGC MDD GWAS summary statistics file
(../inputs/gwas_mdd.tsv) in bounded chunks on local macOS to extract regional SNP windows
around:
  1. The three coloc-ready loci from ../inputs/coloc_ready_loci.tsv
     (locus_0060 DNAH10, locus_0208 SYNE1, locus_0234 RAB2A) +/- 250 kb
  2. Priority SMR / cross-method TWAS causal candidates (DCC, SLC12A5, PPP6C,
     GPX1, SCAMP1, plus HEIDI-failing ANKK1 negative control) using verified GRCh37
     gene coordinates from ../inputs/gene_level_cross_trait_all.tsv.gz +/- 250 kb.

Computes exact Wakefield Asymptotic Bayes Factors (ABF), posterior inclusion
probabilities (PIP), 95%/99% credible sets, prior-variance sensitivity (W = 0.0225,
0.04, 0.09), and LocusZoom-style regional Manhattan + PIP plots.
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


DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent


def current_rss_mb() -> float:
    """Return process resident set size in MiB."""
    return psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)


def wakefield_abf(beta: np.ndarray, se: np.ndarray, W: float = 0.04, snp_prior: float = 1e-4) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute Wakefield log(ABF), normalized single-variant PIPs, and locus PP(H1)."""
    v = np.square(se)
    z = beta / se
    r = W / (W + v)
    log_abf = 0.5 * (np.log(1.0 - r) + r * np.square(z))
    m = float(np.max(log_abf))
    weights = np.exp(log_abf - m)
    pip = weights / np.sum(weights)
    # Locus-level posterior probability of having a causal signal under prior pi_i = snp_prior
    sum_abf_log = m + float(np.log(np.sum(weights)))
    log_odds_h1 = np.log(snp_prior) + sum_abf_log
    pp_h1 = float(1.0 / (1.0 + np.exp(-np.clip(log_odds_h1, -500, 500))))
    return log_abf, pip, pp_h1


def credible_set_size(pips_sorted: np.ndarray, threshold: float) -> int:
    """Return minimum number of ranked SNPs whose cumulative PIP >= threshold."""
    csum = np.cumsum(pips_sorted)
    idx = int(np.searchsorted(csum, threshold, side="left"))
    return min(idx + 1, len(pips_sorted))


def build_target_windows(coloc_loci: pd.DataFrame, magma: pd.DataFrame, flank_bp: int = 250_000) -> list[dict]:
    """Combine coloc_ready_loci and key GRCh37 SMR/TWAS gene intervals into non-overlapping scan targets."""
    targets = []
    for _, row in coloc_loci.iterrows():
        targets.append({
            "locus_tag": f"{row['locus_id']}_{row['gene']}",
            "gene": str(row["gene"]),
            "category": "Coloc_Top3_Locus",
            "chr": int(row["chr"]),
            "win_start": max(1, int(row["start"]) - flank_bp),
            "win_end": int(row["end"]) + flank_bp,
        })

    # Verify GRCh37 gene bounds for top SMR/TWAS candidates directly from gene_level_cross_trait_all.tsv.gz
    mag_idx = magma.dropna(subset=["SYMBOL", "CHR_MDD", "START_MDD", "STOP_MDD"]).drop_duplicates("SYMBOL").set_index("SYMBOL")
    extra_genes = [
        ("SMR_Tier1_Causal", "DCC"),
        ("SMR_Tier1_Causal", "SLC12A5"),
        ("SMR_Tier1_Causal", "PPP6C"),
        ("SMR_Tier1_Causal", "GPX1"),
        ("SMR_Tier1_Causal", "SCAMP1"),
        ("SMR_HEIDI_Fail_Ctrl", "ANKK1"),
    ]
    for cat, g in extra_genes:
        if g in mag_idx.index:
            r = mag_idx.loc[g]
            targets.append({
                "locus_tag": f"GRCh37_{g}",
                "gene": g,
                "category": cat,
                "chr": int(r["CHR_MDD"]),
                "win_start": max(1, int(r["START_MDD"]) - flank_bp),
                "win_end": int(r["STOP_MDD"]) + flank_bp,
            })
    return targets


def stream_extract_gwas_loci(gwas_path: Path, targets: list[dict], chunksize: int = 500_000) -> tuple[dict[str, pd.DataFrame], dict]:
    """Stream 463 MB gwas_mdd.tsv safely and collect SNPs falling within each target window."""
    by_chr: dict[int, list[dict]] = {}
    for t in targets:
        by_chr.setdefault(t["chr"], []).append(t)

    buffers: dict[str, list[pd.DataFrame]] = {t["locus_tag"]: [] for t in targets}
    total_snps = 0
    peak_rss = current_rss_mb()

    t0 = time.perf_counter()
    for chunk in pd.read_csv(gwas_path, sep="\t", chunksize=chunksize):
        total_snps += len(chunk)
        # Filter chunk first to relevant chromosomes
        chrs_present = set(chunk["chr"].unique()) & set(by_chr.keys())
        if not chrs_present:
            continue
        for c in chrs_present:
            c_sub = chunk[chunk["chr"] == c]
            for t in by_chr[c]:
                hit = c_sub[(c_sub["ps"] >= t["win_start"]) & (c_sub["ps"] <= t["win_end"])]
                if len(hit):
                    buffers[t["locus_tag"]].append(hit)
        peak_rss = max(peak_rss, current_rss_mb())
    elapsed = time.perf_counter() - t0

    extracted = {}
    for t in targets:
        tag = t["locus_tag"]
        df = pd.concat(buffers[tag], ignore_index=True) if buffers[tag] else pd.DataFrame()
        if len(df):
            df = df.dropna(subset=["beta", "se", "p_wald"])
            df = df[df["se"] > 0].drop_duplicates(subset=["rs"]).sort_values("ps").reset_index(drop=True)
        extracted[tag] = df

    bench = {
        "gwas_file_resolved": str(gwas_path.resolve()),
        "total_gwas_snps_streamed": int(total_snps),
        "stream_elapsed_seconds": round(elapsed, 3),
        "peak_rss_during_stream_mb": round(peak_rss, 2),
        "chunksize_rows": chunksize,
    }
    return extracted, bench


def plot_regional_finemap(extracted: dict[str, pd.DataFrame], summary_df: pd.DataFrame, fig_dir: Path) -> None:
    """Create a 6-panel LocusZoom + Wakefield PIP figure for top loci."""
    panel_tags = [
        "locus_0208_SYNE1",
        "locus_0060_DNAH10",
        "locus_0234_RAB2A",
        "GRCh37_DCC",
        "GRCh37_SLC12A5",
        "GRCh37_PPP6C",
    ]
    present = [t for t in panel_tags if t in extracted and len(extracted[t]) > 0]
    if not present:
        return

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8), dpi=300)
    axes_flat = axes.flatten()
    sum_idx = summary_df.set_index("locus_tag")

    for ax, tag in zip(axes_flat, present):
        df = extracted[tag].copy()
        df["neglog10p"] = -np.log10(df["p_wald"].clip(lower=1e-300))
        df["pos_mb"] = df["ps"] / 1e6
        in_cs = df["in_95cs"].to_numpy(dtype=bool)

        ax.scatter(df.loc[~in_cs, "pos_mb"], df.loc[~in_cs, "neglog10p"], s=10, color="#7D8B94", alpha=0.55, label="Other SNPs")
        ax.scatter(df.loc[in_cs, "pos_mb"], df.loc[in_cs, "neglog10p"], s=24, color="#C66B45", edgecolor="#233746", linewidth=0.5, zorder=4, label="95% Credible Set")

        top_row = df.loc[df["PIP_W004"].idxmax()]
        ax.scatter([top_row["pos_mb"]], [top_row["neglog10p"]], s=65, marker="D", color="#187C92", edgecolor="#233746", linewidth=0.7, zorder=5)
        ax.annotate(
            f"{top_row['rs']}\nPIP={top_row['PIP_W004']:.2f}",
            xy=(top_row["pos_mb"], top_row["neglog10p"]),
            xytext=(10, 4) if "PPP6C" in tag or "DNAH10" in tag else (8, -4),
            textcoords="offset points",
            fontsize=7.2,
            color="#233746",
            fontweight="bold",
        )
        ax.set_ylim(bottom=0, top=max(10.5, float(df["neglog10p"].max()) * 1.15))
        ax.axhline(-np.log10(5e-8), color="#B83B2A", linestyle="--", linewidth=0.9)

        meta = sum_idx.loc[tag]
        ax.set_title(
            f"{meta['gene']} (chr{int(meta['chr'])})\nCS95={int(meta['cs95_size'])}/{int(meta['n_snps'])} SNPs | minP={meta['min_p']:.2e}",
            fontsize=9.5,
            fontweight="bold",
            loc="left",
            color="#233746",
        )
        ax.set_xlabel("Genomic Position (Mb, GRCh37)")
        ax.set_ylabel(r"$-\log_{10}(P_{\mathrm{MDD}})$")
        if tag == present[0]:
            ax.legend(frameon=False, fontsize=7.5, loc="upper left")

    fig.tight_layout(pad=1.5)
    fig.savefig(fig_dir / "Fig1_Scheme2_Regional_Manhattan_and_PIP.png", dpi=300)
    plt.close(fig)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stream real GWAS summary stats and run Wakefield ABF fine-mapping.")
    parser.add_argument("--in-dir", type=Path, default=DEFAULT_BASE_DIR / "inputs")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_BASE_DIR / "outputs")
    parser.add_argument("--gwas-path", type=Path, default=None)
    parser.add_argument("--chunksize", type=int, default=500_000)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> tuple[pd.DataFrame, dict]:
    args = parse_args(argv)
    in_dir: Path = args.in_dir
    out_dir: Path = args.out_dir
    gwas_path: Path = args.gwas_path if args.gwas_path is not None else (in_dir / "gwas_mdd.tsv")
    if not gwas_path.exists():
        raise FileNotFoundError(f"Missing GWAS summary file: {gwas_path}")

    loc_dir = out_dir / "loci_extracted"
    tab_dir = out_dir / "tables"
    fig_dir = out_dir / "figures"
    bm_dir = out_dir / "benchmarks"
    for d in (loc_dir, tab_dir, fig_dir, bm_dir):
        d.mkdir(parents=True, exist_ok=True)

    coloc_loci = pd.read_csv(in_dir / "coloc_ready_loci.tsv", sep="\t")
    magma = pd.read_csv(in_dir / "gene_level_cross_trait_all.tsv.gz", sep="\t")

    targets = build_target_windows(coloc_loci, magma, flank_bp=250_000)
    extracted, bench = stream_extract_gwas_loci(gwas_path, targets, chunksize=args.chunksize)

    summary_rows = []
    cs_frames = []

    for t in targets:
        tag = t["locus_tag"]
        df = extracted[tag]
        if len(df) == 0:
            continue

        beta = df["beta"].to_numpy(dtype=float)
        se = df["se"].to_numpy(dtype=float)
        _, pip_04, pp_h1_04 = wakefield_abf(beta, se, W=0.04)
        _, pip_0225, _ = wakefield_abf(beta, se, W=0.0225)
        _, pip_09, _ = wakefield_abf(beta, se, W=0.09)

        df["PIP_W004"] = pip_04
        df["PIP_W00225"] = pip_0225
        df["PIP_W009"] = pip_09
        df = df.sort_values(["PIP_W004", "p_wald"], ascending=[False, True]).reset_index(drop=True)
        df["cum_PIP"] = df["PIP_W004"].cumsum()
        cs95_n = credible_set_size(df["PIP_W004"].to_numpy(), 0.95)
        cs99_n = credible_set_size(df["PIP_W004"].to_numpy(), 0.99)
        df["in_95cs"] = np.arange(len(df)) < cs95_n
        df["in_99cs"] = np.arange(len(df)) < cs99_n
        df["locus_tag"] = tag
        df["target_gene"] = t["gene"]

        extracted[tag] = df
        df.sort_values("ps").to_csv(loc_dir / f"{tag}_gwas_grch37.tsv", sep="\t", index=False)

        top = df.iloc[0]
        top_w0225 = df.loc[df["PIP_W00225"].idxmax(), "rs"]
        top_w09 = df.loc[df["PIP_W009"].idxmax(), "rs"]
        cs95_df = df[df["in_95cs"]].copy()
        cs_frames.append(cs95_df)

        summary_rows.append({
            "locus_tag": tag,
            "gene": t["gene"],
            "category": t["category"],
            "chr": t["chr"],
            "win_start": t["win_start"],
            "win_end": t["win_end"],
            "n_snps": len(df),
            "lead_snp": top["rs"],
            "lead_pos": int(top["ps"]),
            "min_p": float(top["p_wald"]),
            "lead_beta": float(top["beta"]),
            "lead_se": float(top["se"]),
            "lead_PIP_W004": float(top["PIP_W004"]),
            "locus_PP_H1": float(pp_h1_04),
            "cs95_size": int(cs95_n),
            "cs99_size": int(cs99_n),
            "prior_robust_top_snp_same": bool(top["rs"] == top_w0225 == top_w09),
        })

    sum_df = pd.DataFrame(summary_rows).sort_values("min_p").reset_index(drop=True)
    cs_all = pd.concat(cs_frames, ignore_index=True) if cs_frames else pd.DataFrame()

    sum_df.to_csv(tab_dir / "finemap_locus_summary.tsv", sep="\t", index=False)
    cs_all.to_csv(tab_dir / "finemap_credible_sets_95pct.tsv", sep="\t", index=False)
    plot_regional_finemap(extracted, sum_df, fig_dir)

    with open(bm_dir / "gwas_stream_finemap_benchmark.json", "w") as f:
        json.dump(bench, f, indent=2)

    print("Scheme 2 Step 1 (Real GWAS Extraction + Wakefield ABF Fine-Mapping) finished.")
    print(json.dumps(bench, indent=2))
    if len(sum_df):
        print(sum_df[["locus_tag", "n_snps", "lead_snp", "min_p", "lead_PIP_W004", "cs95_size", "prior_robust_top_snp_same"]].to_string(index=False))
    return sum_df, bench


if __name__ == "__main__":
    main()
