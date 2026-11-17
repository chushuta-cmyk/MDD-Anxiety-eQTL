#!/usr/bin/env bash
# Real BrainMeta / eQTLGen SMR + HEIDI Execution Script (Scheme 2)
# Requires external Linux/macOS SMR binary + BrainMeta BESD + 1000G EUR LD panel.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
IN_DIR="$BASE_DIR/inputs"
EQTL_DIR="$IN_DIR/external_eqtl"
OUT_DIR="$BASE_DIR/outputs/tables"

GWAS_FILE="$IN_DIR/gwas_mdd.tsv"
GWAS_COJO="$BASE_DIR/outputs/loci_extracted/gwas_mdd.ma"
EQTL_BESD="$EQTL_DIR/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr"
BFILE_1KG="$EQTL_DIR/1000G_EUR_Phase3_GRCh37/EUR"

CHECK_ONLY=0
if [[ "${1:-}" == "--check-only" ]]; then
  CHECK_ONLY=1
fi

for req in "$GWAS_FILE"; do
  [[ -e "$req" ]] || { echo "Missing $req"; exit 1; }
done

HEADER="$(head -n 1 "$GWAS_FILE" | tr '\t' ' ')"
if [[ "$HEADER" != "chr rs ps allele1 allele0 af beta se p_wald" ]]; then
  echo "Unexpected gwas_mdd.tsv header: $HEADER"
  exit 1
fi

HAS_SMR=0
if command -v smr >/dev/null 2>&1; then HAS_SMR=1; fi
HAS_BFILE=0
if [[ -f "${BFILE_1KG}.bed" ]]; then HAS_BFILE=1; fi

if [[ "$CHECK_ONLY" -eq 1 ]]; then
  echo "SMR_PREFLIGHT_OK: gwas_schema_valid=1, smr_cli=${HAS_SMR}, plink_1kg=${HAS_BFILE}"
  exit 0
fi

if [[ "$HAS_SMR" -eq 0 ]]; then
  echo "ERROR: 'smr' executable not found in PATH."
  echo "Download SMR binary from https://yanglab.westlake.edu.cn/software/smr/#Download"
  exit 1
fi

if [[ "$HAS_BFILE" -eq 0 ]]; then
  echo "Missing 1000G EUR PLINK files at ${BFILE_1KG}.{bed,bim,fam}"
  echo "Place downloaded 1000G EUR LD reference & BrainMeta BESD under $EQTL_DIR"
  exit 1
fi

# Convert real gwas_mdd.tsv (chr, rs, ps, allele1, allele0, af, beta, se, p_wald) to COJO .ma format
if [[ ! -f "$GWAS_COJO" ]]; then
  echo "Formatting real gwas_mdd.tsv -> $GWAS_COJO ..."
  awk 'BEGIN{OFS="\t"; print "SNP","A1","A2","freq","b","se","p","N"}
       NR>1 && $2!="" && $2!="." && $7!="NA" && $8+0>0 && $9+0>0 && !seen[$2]++ {print $2,$4,$5,$6,$7,$8,$9,500000}' "$GWAS_FILE" > "$GWAS_COJO"
fi

echo "Running SMR + HEIDI across autosomes 1..22 ..."
for chr in $(seq 1 22); do
  besd_chr="${EQTL_BESD}${chr}"
  if [[ -f "${besd_chr}.besd" ]]; then
    smr \
      --bfile "$BFILE_1KG" \
      --gwas-summary "$GWAS_COJO" \
      --beqtl-summary "$besd_chr" \
      --out "$OUT_DIR/mdd_brainmeta_smr_chr${chr}" \
      --thread-num 4 \
      --peqtl-smr 5e-8 \
      --peqtl-heidi 1.57e-3 \
      --heidi-min-m 5
  fi
done
echo "SMR + HEIDI finished."
