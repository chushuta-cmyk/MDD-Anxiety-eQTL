#!/usr/bin/env python3
"""Audit Full 16-Panel Brain SMR+HEIDI Results Across Prioritized Loci & Mitochondrial Genes (Scheme 2, Step 5).

Reads real genome-wide BrainMeta + GTEx v8/v10 brain eQTL SMR merged outputs staged under
../inputs/smr_16panel_archive/ (MDD 51,994 rows; Anxiety 38,458 rows):
  1. Evaluates cross-panel SMR significance & HEIDI stability across 16 brain eQTL panels
     for prioritized GWAS/SMR candidates (SLC12A5, DCC, PPP6C, GPX1, AMT, SLC25A12, CDK2AP1,
     ANKK1, NEK4, MRM2, SYNE1, RAB2A, DNAH10).
  2. Scans positional genomic neighborhoods around DNAH10, SYNE1, and RAB2A (+/- 1 Mb) to check
     whether alternative nearby genes carry brain eQTL causal evidence (identifying CDK2AP1 near DNAH10).
  3. Audits all 148 GO:0006119 OXPHOS genes, 100 MitoCarta core OXPHOS subunits, and mitochondrial
     carriers (SLC25A* family) across all 16 brain eQTL panels in MDD and Anxiety.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_BASE_DIR = Path(__file__).resolve().parent.parent


def summarize_candidate_replication(mdd_smr: pd.DataFrame, anx_smr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Summarize multi-panel SMR significance & HEIDI >= 0.05 consistency per candidate gene."""
    rows = []
    for g in genes:
        sub_m = mdd_smr[mdd_smr["index"] == g]
        sub_a = anx_smr[anx_smr["index"] == g]

        m_sig = sub_m[sub_m["p_SMR"] < 5e-5]
        a_sig = sub_a[sub_a["p_SMR"] < 5e-5]
        m_pass = m_sig[m_sig["p_HEIDI"] >= 0.05]
        a_pass = a_sig[a_sig["p_HEIDI"] >= 0.05]

        best_m = sub_m.sort_values("p_SMR").iloc[0] if len(sub_m) else None
        best_a = sub_a.sort_values("p_SMR").iloc[0] if len(sub_a) else None
        sign_concordant_mdd = (
            bool((m_sig["b_SMR"] > 0).all() or (m_sig["b_SMR"] < 0).all())
            if len(m_sig) > 0 else np.nan
        )

        rows.append({
            "Gene": g,
            "MDD_total_QTL_panels": int(len(sub_m)),
            "MDD_panels_pSMR_lt_5e5": int(len(m_sig)),
            "MDD_panels_SMR_and_HEIDI_pass": int(len(m_pass)),
            "MDD_sign_concordant_sig": sign_concordant_mdd,
            "MDD_best_QTL": str(best_m["qtl_name"]) if best_m is not None else "NA",
            "MDD_best_pSMR": float(best_m["p_SMR"]) if best_m is not None else np.nan,
            "MDD_best_pHEIDI": float(best_m["p_HEIDI"]) if best_m is not None else np.nan,
            "MDD_best_bSMR": float(best_m["b_SMR"]) if best_m is not None else np.nan,
            "ANX_total_QTL_panels": int(len(sub_a)),
            "ANX_panels_pSMR_lt_5e5": int(len(a_sig)),
            "ANX_panels_SMR_and_HEIDI_pass": int(len(a_pass)),
            "ANX_best_QTL": str(best_a["qtl_name"]) if best_a is not None else "NA",
            "ANX_best_pSMR": float(best_a["p_SMR"]) if best_a is not None else np.nan,
            "ANX_best_pHEIDI": float(best_a["p_HEIDI"]) if best_a is not None else np.nan,
        })
    return pd.DataFrame(rows).sort_values("MDD_best_pSMR", na_position="last").reset_index(drop=True)


def scan_positional_windows(mdd_smr: pd.DataFrame) -> pd.DataFrame:
    """Scan +/- 1 Mb neighborhoods around Coloc Top3 positional genes for alternative causal SMR genes."""
    windows = [
        ("locus_0208_SYNE1_chr6", 6, 151_200_000, 153_200_000),
        ("locus_0234_RAB2A_chr8", 8, 59_500_000, 61_500_000),
        ("locus_0060_DNAH10_chr12", 12, 122_600_000, 124_800_000),
    ]
    frames = []
    for tag, chrom, start, end in windows:
        sub = mdd_smr[
            (mdd_smr["chr"] == chrom)
            & (mdd_smr["Probe_bp"] >= start)
            & (mdd_smr["Probe_bp"] <= end)
            & (mdd_smr["p_SMR"] < 1e-3)
        ].sort_values("p_SMR").copy()
        if len(sub):
            sub.insert(0, "window_tag", tag)
            frames.append(sub[["window_tag", "index", "qtl_name", "ProbeChr", "Probe_bp", "topSNP", "p_SMR", "p_HEIDI", "b_SMR"]])
    if frames:
        return pd.concat(frames, ignore_index=True)
    return pd.DataFrame(columns=["window_tag", "index", "qtl_name", "ProbeChr", "Probe_bp", "topSNP", "p_SMR", "p_HEIDI", "b_SMR"])


def audit_mitochondrial_smr(
    mdd_smr: pd.DataFrame,
    anx_smr: pd.DataFrame,
    ox_genes: set[str],
    mito_genes: set[str],
) -> pd.DataFrame:
    """Check GO:0006119 (148), MitoCarta_core (100), and SLC25A* carriers across 16 brain eQTL panels."""
    out_frames = []
    for trait_name, df in [("MDD", mdd_smr), ("ANX", anx_smr)]:
        is_ox = df["index"].isin(ox_genes)
        is_mito = df["index"].isin(mito_genes)
        is_slc25 = df["index"].astype(str).str.startswith("SLC25A")
        sub = df[(is_ox | is_mito | is_slc25) & (df["p_SMR"] < 1e-3)].copy()
        if len(sub):
            sub.insert(0, "trait", trait_name)
            sub["in_GO0006119"] = sub["index"].isin(ox_genes)
            sub["in_MitoCarta_core100"] = sub["index"].isin(mito_genes)
            sub["is_SLC25A_carrier"] = sub["index"].astype(str).str.startswith("SLC25A")
            out_frames.append(
                sub[[
                    "trait", "index", "in_GO0006119", "in_MitoCarta_core100", "is_SLC25A_carrier",
                    "qtl_name", "topSNP", "b_SMR", "p_SMR", "p_HEIDI", "GWAS_LOCUS",
                ]].sort_values("p_SMR")
            )
    if out_frames:
        return pd.concat(out_frames, ignore_index=True)
    return pd.DataFrame()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit 16-panel BrainMeta+GTEx SMR archive across candidates and mitochondrial genes.")
    parser.add_argument("--in-dir", type=Path, default=DEFAULT_BASE_DIR / "inputs")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_BASE_DIR / "outputs")
    parser.add_argument("--mdd-esmr", type=Path, default=None)
    parser.add_argument("--anx-esmr", type=Path, default=None)
    parser.add_argument("--allow-missing", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> dict:
    args = parse_args(argv)
    in_dir: Path = args.in_dir
    out_dir: Path = args.out_dir
    tab_dir = out_dir / "tables"
    bm_dir = out_dir / "benchmarks"
    for d in (tab_dir, bm_dir):
        d.mkdir(parents=True, exist_ok=True)

    mdd_path = args.mdd_esmr if args.mdd_esmr else (in_dir / "smr_16panel_archive" / "MDD_trait_eSMR.merged.tsv")
    anx_path = args.anx_esmr if args.anx_esmr else (in_dir / "smr_16panel_archive" / "ANX_trait_eSMR.merged.tsv")

    if not (mdd_path.exists() and anx_path.exists()):
        msg = f"16-panel SMR archive not found ({mdd_path}, {anx_path})."
        if args.allow_missing:
            print("SKIP:", msg)
            return {"status": "skipped_missing_archive"}
        raise FileNotFoundError(msg)

    t0 = time.perf_counter()
    mdd_smr = pd.read_csv(mdd_path, sep="\t")
    anx_smr = pd.read_csv(anx_path, sep="\t")

    ox_path = DEFAULT_BASE_DIR.parent / "01_Scheme1_snRNA_Spatial_scDRS" / "inputs" / "OXPHOS_GO0006119_148genes.tsv"
    if not ox_path.exists():
        ox_path = in_dir / "OXPHOS_GO0006119_148genes.tsv"
    ox_genes = set(pd.read_csv(ox_path, sep="\t")["Gene"].astype(str).str.strip()) if ox_path.exists() else set()

    mito_path = in_dir / "MitoCarta_core_frozen.tsv"
    mito_genes = (
        set(pd.read_csv(mito_path, sep="\t")["MitoCartaSymbol"].astype(str).str.strip().str.replace('"', "", regex=False))
        if mito_path.exists()
        else set()
    )

    genes = ["SLC12A5", "DCC", "PPP6C", "GPX1", "AMT", "SLC25A12", "CDK2AP1", "ANKK1", "NEK4", "MRM2", "SYNE1", "RAB2A", "DNAH10"]
    rep_df = summarize_candidate_replication(mdd_smr, anx_smr, genes)
    win_df = scan_positional_windows(mdd_smr)
    mito_df = audit_mitochondrial_smr(mdd_smr, anx_smr, ox_genes, mito_genes)

    rep_df.to_csv(tab_dir / "smr_16panel_prioritized_replication.tsv", sep="\t", index=False)
    win_df.to_csv(tab_dir / "smr_16panel_positional_neighbor_scan.tsv", sep="\t", index=False)
    mito_df.to_csv(tab_dir / "smr_16panel_oxphos_mitocarta_audit.tsv", sep="\t", index=False)

    ox_mdd_strict = int(((mdd_smr["index"].isin(ox_genes)) & (mdd_smr["p_SMR"] < 1e-4)).sum())
    summary = {
        "script": "05_audit_16panel_brain_smr_and_oxphos.py",
        "elapsed_seconds": round(time.perf_counter() - t0, 3),
        "mdd_rows": int(len(mdd_smr)),
        "anx_rows": int(len(anx_smr)),
        "mdd_qtl_panels": int(mdd_smr["qtl_name"].nunique()),
        "go0006119_mdd_smr_lt_1e4": ox_mdd_strict,
        "mito_hits_lt_1e3": int(len(mito_df)),
        "neighbor_hits_lt_1e3": int(len(win_df)),
    }
    with open(bm_dir / "smr_16panel_audit_benchmark.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Scheme 2 Step 5 (16-Panel BrainMeta+GTEx SMR Audit) finished.")
    print(json.dumps(summary, indent=2))
    print(rep_df[["Gene", "MDD_panels_pSMR_lt_5e5", "MDD_panels_SMR_and_HEIDI_pass", "MDD_best_pSMR", "MDD_best_pHEIDI"]].to_string(index=False))
    return summary


if __name__ == "__main__":
    main()
