#!/usr/bin/env python3
"""Harmonize Bryois 2022 cell-type eQTL effect alleles against MDD & ANX GWAS risk alleles."""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Allele harmonization for Bryois eQTL colocalization")
    p.add_argument("--snp-pos", type=Path, default=SCRIPT_DIR / "external_downloads" / "bryois_snp_pos.txt.gz")
    p.add_argument("--coloc-tsv", type=Path, default=SCRIPT_DIR / "outputs" / "tables" / "aim3_bryois_8celltype_coloc_MDD_and_ANX.tsv")
    p.add_argument("--loci-dir", type=Path, default=SCRIPT_DIR / "outputs" / "loci_dualtrait")
    p.add_argument("--out-tsv", type=Path, default=SCRIPT_DIR / "outputs" / "tables" / "p0_bryois_allelic_harmonized_all.tsv")
    p.add_argument("--check-only", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.check_only:
        print(f"P0_ALLELE_PREFLIGHT: coloc={int(args.coloc_tsv.exists())}, loci={len(list(args.loci_dir.glob('*_dualtrait_grch37.tsv')))}")
        return 0
    bdf = pd.read_csv(args.coloc_tsv, sep="\t")
    need_snps = set(bdf["top_eqtl_snp"].dropna().astype(str))
    snp_alleles = {}
    with gzip.open(args.snp_pos, "rt") as f:
        hdr = f.readline().strip().split("\t")
        ci = {k: i for i, k in enumerate(hdr)}
        for line in f:
            toks = line.rstrip("\n").split("\t")
            rs = toks[ci["SNP"]]
            if rs in need_snps:
                snp_alleles[rs] = (toks[ci["effect_allele"]].upper(), toks[ci["other_allele"]].upper())
    gw_map = {}
    for lf in args.loci_dir.glob("*_dualtrait_grch37.tsv"):
        df = pd.read_csv(lf, sep="\t")
        for _, r in df[df["SNP"].isin(need_snps)].iterrows():
            gw_map[r["SNP"]] = {
                "A1_MDD": str(r["A1_MDD"]).upper(), "A2_MDD": str(r["A2_MDD"]).upper(),
                "BETA_MDD": float(r["BETA_MDD"]), "P_MDD": float(r["P_MDD"]),
                "A1_ANX": str(r["A1_ANX"]).upper(), "A2_ANX": str(r["A2_ANX"]).upper(),
                "BETA_ANX": float(r["BETA_ANX"]), "P_ANX": float(r["P_ANX"]),
            }
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rows = []
    for _, r in bdf.iterrows():
        rs = str(r["top_eqtl_snp"])
        if rs not in snp_alleles or rs not in gw_map:
            continue
        ea, oa = snp_alleles[rs]
        gw = gw_map[rs]
        tr = r["trait"]
        a1, a2, bg, pg = gw[f"A1_{tr}"], gw[f"A2_{tr}"], gw[f"BETA_{tr}"], gw[f"P_{tr}"]
        if (ea, oa) == (a1, a2) or (comp.get(ea), comp.get(oa)) == (a1, a2):
            be, status = float(r["top_eqtl_beta"]), "aligned"
        elif (ea, oa) == (a2, a1) or (comp.get(ea), comp.get(oa)) == (a2, a1):
            be, status = -float(r["top_eqtl_beta"]), "flipped"
        else:
            be, status = np.nan, "mismatch"
        risk_a, oth_a = (a2, a1) if bg < 0 else (a1, a2)
        eq_risk = -be if bg < 0 else be
        gw_risk = abs(bg)
        rows.append({
            "locus": r["locus"], "cell_type": r["cell_type"], "gene": r["gene"], "trait": tr,
            "PP.H4": r["PP.H4"], "top_eqtl_snp": rs, "risk_allele": risk_a, "other_allele": oth_a,
            "align_status": status, "gwas_beta_risk": gw_risk, "gwas_p": pg,
            "eqtl_beta_on_risk_allele": eq_risk, "top_eqtl_p": r["top_eqtl_p"],
            "risk_allele_effect_on_expr": "INCREASES_EXPR" if eq_risk > 0 else ("DECREASES_EXPR" if eq_risk < 0 else "NA"),
        })
    pd.DataFrame(rows).to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[SUCCESS] Harmonized {len(rows)} Bryois coloc rows to risk-increasing alleles.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
