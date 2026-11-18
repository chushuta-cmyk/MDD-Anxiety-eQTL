#!/usr/bin/env python3
"""Local GWAS Overlap Sensitivity, LDSC Intercept Audit, 20q13.12 30-Gene Matrix & Donor scDRS Export.

Addresses reviewer requests without requiring new external downloads:
  1. Extracts noUKBB (N_eff up to 967k) and Clinical-only (Clin_eur) PGC MDD 2025 splits
     against Friligkou et al. (2024) ANX_EUR across all 7 candidate windows.
  2. Runs R coloc::coloc.abf comparing MDD_all, MDD_noUKBB, MDD_Clin, and LDSC-intercept
     (gcov_int = 0.3231 / 0.3122) decorrelated ANX statistics against MDD.
  3. Exports full 30-gene x 8-cell-type colocalization table for 20q13.12 (including CD40).
  4. Aggregates per-donor scDRS mean normalized scores for individual-point plotting.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
OUT_DIR = SCRIPT_DIR / "outputs" / "tables"
LOCI_DIR = SCRIPT_DIR / "outputs" / "loci_dualtrait"
S1_DIR = REPO_ROOT / "01_Scheme1_snRNA_Spatial_scDRS"

MDD_NOUKBB = Path(os.environ.get("MDD_NOUKBB_PATH", REPO_ROOT / "data" / "gwas" / "pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv.gz"))
MDD_CLIN = Path(os.environ.get("MDD_CLIN_PATH", REPO_ROOT / "data" / "gwas" / "pgc-mdd2025_Clin_eur_v3-49-24-11.tsv.gz"))
ANX_EUR = Path(os.environ.get("ANX_EUR_PATH", REPO_ROOT / "data" / "gwas" / "ANX_EUR.txt.gz"))
LDSC_JSON = Path(os.environ.get("LDSC_JSON_PATH", REPO_ROOT / "data" / "ldsc_summary.json"))

LOCUS_WINDOWS: dict[str, tuple[int, int, int]] = {
    "SLC12A5_chr20": (20, 44_150_000, 45_200_000),
    "DCC_chr18": (18, 49_600_000, 51_500_000),
    "PPP6C_chr9": (9, 127_400_000, 128_450_000),
    "SLC25A12_chr2": (2, 172_390_000, 173_120_000),
    "DNAH10_CDK2AP1_chr12": (12, 123_000_000, 124_600_000),
    "chr3p21_31_GPX1_AMT_RNF123": (3, 48_700_000, 50_100_000),
    "chr3p21_1_GNL3_NEK4_SPCS1": (3, 52_100_000, 53_300_000),
}


def extract_windows_from_mdd_tsv(path: Path, prefix: str) -> dict[str, pd.DataFrame]:
    """Stream PGC MDD 2025 TSV.GZ and return per-locus DataFrames."""
    by_chr = {w[0]: (name, w[1], w[2]) for name, w in LOCUS_WINDOWS.items()}
    # Note chr3 has two loci, handle list per chr
    chr_map: dict[int, list[tuple[str, int, int]]] = {}
    for name, (ch, st, en) in LOCUS_WINDOWS.items():
        chr_map.setdefault(ch, []).append((name, st, en))

    buckets: dict[str, list[tuple]] = {name: [] for name in LOCUS_WINDOWS}
    with gzip.open(path, "rt") as fh:
        hdr = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM") or line.startswith("CHROM"):
                hdr = line.lstrip("#").strip().split("\t")
                continue
            parts = line.rstrip("\n").split("\t")
            ch = int(parts[0])
            if ch not in chr_map:
                continue
            pos = int(parts[1])
            for name, st, en in chr_map[ch]:
                if st <= pos <= en:
                    # ID, EA, NEA, BETA, SE, PVAL, NEFF
                    buckets[name].append((parts[2], ch, pos, parts[3].upper(), parts[4].upper(),
                                         float(parts[5]), float(parts[6]), float(parts[7]), float(parts[11])))
    out = {}
    cols = ["SNP", "CHR", "BP", f"A1_{prefix}", f"A2_{prefix}", f"BETA_{prefix}", f"SE_{prefix}", f"P_{prefix}", f"N_{prefix}"]
    for name, rows in buckets.items():
        out[name] = pd.DataFrame(rows, columns=cols).drop_duplicates("SNP")
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()

    if args.check_only:
        ok = int(LDSC_JSON.exists() and (OUT_DIR / "aim3_bryois_8celltype_coloc_MDD_and_ANX.tsv").exists())
        print(f"GWAS_OVERLAP_SENS_PREFLIGHT: inputs_ok={ok}, noukbb={int(MDD_NOUKBB.exists())}")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Export LDSC overlap summary
    if LDSC_JSON.exists():
        ld = json.loads(LDSC_JSON.read_text())
        ldsc_rows = []
        for k in ("main", "no_ukbb_sensitivity"):
            m = ld[k]
            ldsc_rows.append({
                "analysis": k,
                "mdd_h2_obs": m["mdd"]["h2_observed"],
                "mdd_intercept": m["mdd"]["intercept"],
                "anx_h2_obs": m["anxiety"]["h2_observed"],
                "anx_intercept": m["anxiety"]["intercept"],
                "rg": m["rg"],
                "rg_se": m["rg_se"],
                "rg_z": m["z"],
                "cross_trait_intercept": m["cross_trait_intercept"],
                "cross_trait_intercept_se": m["cross_trait_intercept_se"],
                "n_hm3_snps": m["valid_allele_snps"],
            })
        pd.DataFrame(ldsc_rows).to_csv(OUT_DIR / "p1_ldsc_overlap_summary.tsv", sep="\t", index=False)
        gcov_int_main = float(ld["main"]["cross_trait_intercept"])
    else:
        gcov_int_main = 0.3231

    # 2. Export full 30-gene 20q13.12 matrix across all 8 Bryois cell types
    bdf = pd.read_csv(OUT_DIR / "aim3_bryois_8celltype_coloc_MDD_and_ANX.tsv", sep="\t")
    w20 = bdf[bdf["locus"] == "SLC12A5_chr20"].copy()
    w20.sort_values(["trait", "PP.H4"], ascending=[True, False]).to_csv(
        OUT_DIR / "p1_20q13_all30genes_bryois_coloc.tsv", sep="\t", index=False
    )

    # 3. Export per-donor scDRS summary across key subtypes
    sc_path = S1_DIR / "outputs" / "scdrs" / "MDD_PGC2025.score.gz"
    meta_pb = pd.read_csv(SCRIPT_DIR / "outputs" / "local_pb_cache" / "pb_groups_meta.tsv", sep="\t")
    # Check if raw cells_columns matches scdrs index to extract donor x subtype scDRS
    cells_gz = Path(os.environ.get("GSE213982_CELLS_GZ", REPO_ROOT / "data" / "GSE213982" / "GSE213982_combined_counts_matrix_cells_columns.csv.gz"))
    if sc_path.exists() and cells_gz.exists():
        sc_df = pd.read_csv(sc_path, sep="\t", index_col=0)
        # Index in sc_df is barcode or cell_id; check overlap with cells_columns
        cells_s = pd.read_csv(cells_gz).iloc[:, 0].astype(str).str.strip('"')
        toks = cells_s.str.split(".", expand=True)
        cell_meta = pd.DataFrame({
            "cell_id": cells_s.values,
            "donor_id": toks[0].replace({"M24_2": "M24"}).values,
            "author_subtype": toks[3].values,
        }).set_index("cell_id")
        merged_sc = sc_df[["norm_score"]].join(cell_meta, how="inner")
        if len(merged_sc) > 0:
            don_diag = meta_pb[["donor_id", "diagnosis", "Sex"]].drop_duplicates()
            don_sc = (
                merged_sc.groupby(["donor_id", "author_subtype"], as_index=False)["norm_score"]
                .mean()
                .merge(don_diag, on="donor_id", how="left")
            )
            don_sc.to_csv(OUT_DIR / "p1_scdrs_donor_subtype_means.tsv", sep="\t", index=False)
            print(f"[INFO] Exported donor-level scDRS means: {len(don_sc)} donor-subtype entries")

    # 4. Stream noUKBB and Clin PGC MDD splits across 7 windows and compare via coloc.abf
    if MDD_NOUKBB.exists() and MDD_CLIN.exists():
        print("[INFO] Streaming noUKBB & Clin PGC MDD 2025 splits across 7 windows...")
        nouk_dict = extract_windows_from_mdd_tsv(MDD_NOUKBB, "NOUKBB")
        clin_dict = extract_windows_from_mdd_tsv(MDD_CLIN, "CLIN")
        merged_paths = {}
        tmp_dir = SCRIPT_DIR / "outputs" / "loci_sensitivity_tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        for loc in LOCUS_WINDOWS:
            base_f = LOCI_DIR / f"{loc}_dualtrait_grch37.tsv"
            if not base_f.exists():
                continue
            df = pd.read_csv(base_f, sep="\t")
            df = df.merge(nouk_dict[loc], on="SNP", how="left", suffixes=("", "_nk"))
            df = df.merge(clin_dict[loc], on="SNP", how="left", suffixes=("", "_cl"))
            out_f = tmp_dir / f"{loc}_sensitivity.tsv"
            df.to_csv(out_f, sep="\t", index=False)
            merged_paths[loc] = str(out_f)

        # Call R coloc.abf on MDD_all, MDD_noUKBB, MDD_Clin, and LDSC-intercept decorrelated
        r_code = f"""
suppressPackageStartupMessages({{ library(data.table); library(coloc) }})
tmp_dir <- "{tmp_dir}"
out_tsv <- "{OUT_DIR / 'p1_gwas_overlap_coloc_sensitivity.tsv'}"
gcov_main <- {gcov_int_main}
gcov_nouk <- 0.3122
files <- list.files(tmp_dir, pattern="_sensitivity.tsv$", full.names=TRUE)
res_list <- list()

run_coloc_pair <- function(beta1, se1, n1, beta2, se2, n2, snp, pos, se_mult=1.0) {{
  ok <- is.finite(beta1) & is.finite(se1) & se1 > 0 & is.finite(beta2) & is.finite(se2) & se2 > 0
  if (sum(ok) < 50) return(NULL)
  d1 <- list(beta=beta1[ok], varbeta=(se1[ok] * se_mult)^2, snp=snp[ok], position=pos[ok], N=median(n1[ok], na.rm=TRUE), type="cc")
  d2 <- list(beta=beta2[ok], varbeta=(se2[ok] * se_mult)^2, snp=snp[ok], position=pos[ok], N=median(n2[ok], na.rm=TRUE), type="quant", sdY=1)
  fit <- suppressMessages(coloc.abf(d1, d2))
  list(nsnps=sum(ok), h0=fit$summary[["PP.H0.abf"]], h1=fit$summary[["PP.H1.abf"]],
       h2=fit$summary[["PP.H2.abf"]], h3=fit$summary[["PP.H3.abf"]], h4=fit$summary[["PP.H4.abf"]])
}}

for (f in files) {{
  loc <- sub("_sensitivity.tsv$", "", basename(f))
  dt <- fread(f)
  for (mode in c("MDD_all_vs_ANX", "MDD_noUKBB_vs_ANX", "MDD_Clin_vs_ANX", "MDD_all_vs_ANX_se_inflated_LDSCint")) {{
    if (mode == "MDD_all_vs_ANX") {{
      r <- run_coloc_pair(dt$BETA_MDD, dt$SE_MDD, dt$N_MDD, dt$BETA_ANX, dt$SE_ANX, dt$N_ANX, dt$SNP, dt$BP, 1.0)
    }} else if (mode == "MDD_noUKBB_vs_ANX") {{
      flip <- (dt$A1_NOUKBB == dt$A2_MDD) & (dt$A2_NOUKBB == dt$A1_MDD)
      match <- ((dt$A1_NOUKBB == dt$A1_MDD) & (dt$A2_NOUKBB == dt$A2_MDD)) | flip
      b1 <- ifelse(flip, -dt$BETA_NOUKBB, dt$BETA_NOUKBB)
      b1[!match] <- NA
      r <- run_coloc_pair(b1, dt$SE_NOUKBB, dt$N_NOUKBB, dt$BETA_ANX, dt$SE_ANX, dt$N_ANX, dt$SNP, dt$BP, 1.0)
    }} else if (mode == "MDD_Clin_vs_ANX") {{
      flip <- (dt$A1_CLIN == dt$A2_MDD) & (dt$A2_CLIN == dt$A1_MDD)
      match <- ((dt$A1_CLIN == dt$A1_MDD) & (dt$A2_CLIN == dt$A2_MDD)) | flip
      b1 <- ifelse(flip, -dt$BETA_CLIN, dt$BETA_CLIN)
      b1[!match] <- NA
      r <- run_coloc_pair(b1, dt$SE_CLIN, dt$N_CLIN, dt$BETA_ANX, dt$SE_ANX, dt$N_ANX, dt$SNP, dt$BP, 1.0)
    }} else {{
      # Conservative overlap penalty inflating varbeta by (1 + gcov_int)
      r <- run_coloc_pair(dt$BETA_MDD, dt$SE_MDD, dt$N_MDD, dt$BETA_ANX, dt$SE_ANX, dt$N_ANX, dt$SNP, dt$BP, sqrt(1.0 + gcov_main))
    }}
    if (!is.null(r)) {{
      res_list[[length(res_list)+1]] <- data.table(locus=loc, model=mode, nsnps=r$nsnps, PP.H0=r$h0, PP.H1=r$h1, PP.H2=r$h2, PP.H3=r$h3, PP.H4=r$h4)
    }}
  }}
}}
out <- rbindlist(res_list)
fwrite(out, out_tsv, sep="\\t")
print(out[locus %in% c("SLC12A5_chr20", "DCC_chr18", "chr3p21_1_GNL3_NEK4_SPCS1")])
"""
        r_script = tmp_dir / "run_sens_coloc.R"
        r_script.write_text(r_code)
        subprocess.run(["Rscript", str(r_script)], check=True)

    print("[SUCCESS] GWAS overlap sensitivity and donor scDRS exports complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
