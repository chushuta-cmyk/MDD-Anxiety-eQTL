#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Server / Local Pipeline: Fix ANX A1 Frequency via 1000G EUR & Run SMR-multi
# ==============================================================================
# Addresses the exact blocker documented in `ANX_SMR_upload_blocker.txt`:
#   ANX summary stats MAF is minor-allele frequency, not guaranteed A1 frequency.
#   Step 1 computes reference A1/A2 frequencies from 1000G EUR PLINK bfiles,
#          matches & flips freq to match ANX_EUR_GWAShug_aligned3.tsv.gz A1,
#          and writes valid COJO .ma format (SNP A1 A2 freq b se p n).
#   Step 2 Runs multi-SNP SMR + HEIDI (`--smr-multi --heidi-mtd 1`) across
#          BrainMeta cis-eQTL BESD, BrainMeta cis-sQTL BESD, and Bryois BESDs.
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
S2_ARCH="${REPO_ROOT}/02_Scheme2_GWAS_Causal_Coloc_SMR/inputs/large_n_gwas_archive"

ANX_GWAS="${S2_ARCH}/ANX_EUR_GWAShug_aligned3.tsv.gz"
MDD_MA="${S2_ARCH}/MDD_no23andMe_EUR_SMR.ma.gz"
REF_1KG="${SCRIPT_DIR}/external_downloads/ld_1kg_eur/g1000_eur"
BESD_DIR="${SCRIPT_DIR}/external_downloads/brainmeta_besd"
OUT_DIR="${SCRIPT_DIR}/outputs/smr_multi"
CHECK_ONLY=0
THREADS="${SLURM_CPUS_PER_TASK:-4}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --check-only) CHECK_ONLY=1; shift ;;
    --ref-1kg) REF_1KG="$2"; shift 2 ;;
    --besd-dir) BESD_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "${OUT_DIR}"

HAS_PLINK=0; (command -v plink >/dev/null 2>&1 || command -v plink2 >/dev/null 2>&1) && HAS_PLINK=1
HAS_SMR=0; (command -v smr >/dev/null 2>&1 || command -v smr-1.3.1 >/dev/null 2>&1) && HAS_SMR=1
HAS_ANX=0; [[ -s "${ANX_GWAS}" ]] && HAS_ANX=1
HAS_MDD=0; [[ -s "${MDD_MA}" ]] && HAS_MDD=1
HAS_1KG=0; [[ -f "${REF_1KG}.bed" ]] && HAS_1KG=1

echo "AIM2_SMR_PREFLIGHT: anx_aligned3=${HAS_ANX}, mdd_ma=${HAS_MDD}, plink=${HAS_PLINK}, smr=${HAS_SMR}, bfile_1kg=${HAS_1KG}"

if [[ "${CHECK_ONLY}" -eq 1 ]]; then
  exit 0
fi

if [[ "${HAS_PLINK}" -eq 0 || "${HAS_SMR}" -eq 0 || "${HAS_1KG}" -eq 0 ]]; then
  echo "ERROR: Require plink/plink2, smr, and ${REF_1KG}.{bed,bim,fam} before full run." >&2
  exit 2
fi

SMR_BIN="$(command -v smr || command -v smr-1.3.1)"

# Step 1: Compute 1000G EUR reference allele frequencies
if [[ ! -f "${OUT_DIR}/g1000_eur.frq" ]]; then
  plink --bfile "${REF_1KG}" --freq --out "${OUT_DIR}/g1000_eur" --threads "${THREADS}"
fi

# Step 2: Build ANX COJO .ma with verified 1000G EUR A1 frequencies
python3 - <<PY
import gzip, pandas as pd, numpy as np
from pathlib import Path

frq = pd.read_csv("${OUT_DIR}/g1000_eur.frq", sep=r"\s+", usecols=["SNP","A1","A2","MAF"])
anx = pd.read_csv("${ANX_GWAS}", sep=r"\s+", usecols=["SNP","A1","A2","BETA","SE","P","N"])
m = anx.merge(frq, on="SNP", suffixes=("_gwas", "_ref")).drop_duplicates("SNP")
a1g = m["A1_gwas"].astype(str).str.upper()
a2g = m["A2_gwas"].astype(str).str.upper()
a1r = m["A1_ref"].astype(str).str.upper()
a2r = m["A2_ref"].astype(str).str.upper()

same = (a1g == a1r) & (a2g == a2r)
flip = (a1g == a2r) & (a2g == a1r)
m = m[same | flip].copy()
m["freq"] = np.where((m["A1_gwas"].str.upper() == m["A1_ref"].str.upper()), m["MAF"], 1.0 - m["MAF"])
# Exclude extreme palindromic SNPs near 0.5
pal = (a1g + a2g).isin(["AT","TA","CG","GC"]) & (m["freq"].sub(0.5).abs() < 0.08)
m = m[~pal[m.index]].copy()
out = m[["SNP","A1_gwas","A2_gwas","freq","BETA","SE","P","N"]].rename(
    columns={"A1_gwas":"A1","A2_gwas":"A2","BETA":"b","SE":"se","P":"p","N":"n"}
)
out.to_csv("${OUT_DIR}/ANX_EUR_SMR_1000Gfreq.ma", sep="\t", index=False)
print("Wrote aligned ANX .ma SNPs:", len(out))
PY

# Unpack MDD .ma if zipped
MDD_UNZIP="${OUT_DIR}/MDD_no23andMe_EUR_SMR.ma"
if [[ ! -f "${MDD_UNZIP}" ]]; then
  gzip -dc "${MDD_MA}" > "${MDD_UNZIP}"
fi

# Step 3: Execute --smr-multi --heidi-mtd 1 across all BESDs found in BESD_DIR
for besd_file in "${BESD_DIR}"/*.besd; do
  [[ -e "${besd_file}" ]] || continue
  prefix="${besd_file%.besd}"
  bname="$(basename "${prefix}")"
  echo "Running --smr-multi on BESD=${bname} ..."
  "${SMR_BIN}" \
    --bfile "${REF_1KG}" \
    --gwas-summary "${MDD_UNZIP}" \
    --beqtl-summary "${prefix}" \
    --smr-multi \
    --heidi-mtd 1 \
    --peqtl-smr 5e-8 \
    --maf 0.01 \
    --thread-num "${THREADS}" \
    --out "${OUT_DIR}/MDD_${bname}_smr_multi"

  "${SMR_BIN}" \
    --bfile "${REF_1KG}" \
    --gwas-summary "${OUT_DIR}/ANX_EUR_SMR_1000Gfreq.ma" \
    --beqtl-summary "${prefix}" \
    --smr-multi \
    --heidi-mtd 1 \
    --peqtl-smr 5e-8 \
    --maf 0.01 \
    --thread-num "${THREADS}" \
    --out "${OUT_DIR}/ANX_${bname}_smr_multi"
done
