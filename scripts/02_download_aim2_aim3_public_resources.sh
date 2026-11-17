#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Aim 2 & Aim 3 Public Resource Downloader & Disk Safety Guardrail
# ==============================================================================
# Usage:
#   bash 02_download_aim2_aim3_public_resources.sh --check-only   # Check disk & print manifest
#   bash 02_download_aim2_aim3_public_resources.sh --tier p0       # Download GEO/Bryois/SMR binaries
#
# Notes:
#   - Uses official Zenodo record 10.5281/zenodo.5543734 (Bryois 2022 8 brain cell-type cis-eQTLs)
#   - Uses official NCBI GEO FTP for GSE213982 Series Matrix clinical metadata
#   - Prints exact manual/wget instructions for Yang Lab BrainMeta BESD, MetaBrain, and Synapse gated files
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_RAW="${SCRIPT_DIR}/external_downloads"
CHECK_ONLY=0
TIER="p0"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --check-only) CHECK_ONLY=1; shift ;;
    --tier) TIER="${2:-p0}"; shift 2 ;;
    --out-dir) OUT_RAW="$2"; shift 2 ;;
    *) echo "Unknown flag: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "${OUT_RAW}/geo" "${OUT_RAW}/bryois_zenodo_5543734" "${OUT_RAW}/ld_1kg_eur" "${OUT_RAW}/brainmeta_besd" "${OUT_RAW}/bins"

AVAIL_KB="$(df -Pk "${SCRIPT_DIR}" | awk 'NR==2 {print $4}')"
AVAIL_GB=$(( AVAIL_KB / 1024 / 1024 ))
HAS_CURL=0; command -v curl >/dev/null 2>&1 && HAS_CURL=1
HAS_WGET=0; command -v wget >/dev/null 2>&1 && HAS_WGET=1

echo "DOWNLOAD_PREFLIGHT: tier=${TIER}, avail_disk_gb=${AVAIL_GB}, curl=${HAS_CURL}, wget=${HAS_WGET}, target=${OUT_RAW}"

if [[ "${CHECK_ONLY}" -eq 1 ]]; then
  echo "Manifest:"
  echo "  [P0-1] NCBI GEO GSE213982 Series Matrix -> ${OUT_RAW}/geo/"
  echo "         URL: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213982/matrix/"
  echo "  [P0-2] Bryois et al. 2022 8-Cell-Type Brain cis-eQTLs (DOI: 10.5281/zenodo.5543734)"
  echo "         URL: https://zenodo.org/records/5543734"
  echo "  [P0-3] Yang Lab SMR binary (Linux/macOS build from https://github.com/JianYang-Lab/SMR)"
  echo "  [P0-4] Yang Lab BrainMeta v2 cis-eQTL & cis-sQTL BESD (https://yanglab.westlake.edu.cn/software/smr/#DataResource)"
  echo "  [P0-5] 1000 Genomes Phase 3 EUR PLINK bfile (GRCh37) -> ${OUT_RAW}/ld_1kg_eur/"
  echo "  [P1-1] MetaBrain Cortex cis-eQTL (https://www.metabrain.nl/)"
  echo "  [P1-2] ROSMAP Fujita et al. 2024 & Wingo 2021 DLPFC pQTL (Synapse auth required)"
  exit 0
fi

if [[ "${AVAIL_GB}" -lt 20 ]]; then
  echo "ERROR: Free disk space (${AVAIL_GB} GB) < 20 GB minimum safe threshold for uncompressed BESD/eQTL work." >&2
  exit 2
fi

if [[ "${HAS_CURL}" -eq 0 && "${HAS_WGET}" -eq 0 ]]; then
  echo "ERROR: Neither curl nor wget found." >&2
  exit 2
fi

fetch() {
  local url="$1"
  local dest="$2"
  echo "Downloading: ${url} -> ${dest}"
  if [[ "${HAS_CURL}" -eq 1 ]]; then
    curl -fL --retry 3 -o "${dest}" "${url}"
  else
    wget -c -O "${dest}" "${url}"
  fi
}

if [[ "${TIER}" == "p0" || "${TIER}" == "all" ]]; then
  # 1. NCBI GEO GSE213982 Series Matrix (try standard GPL24676 or directory listing first)
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213982/matrix/GSE213982_series_matrix.txt.gz" "${OUT_RAW}/geo/GSE213982_series_matrix.txt.gz" || \
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213982/matrix/GSE213982-GPL24676_series_matrix.txt.gz" "${OUT_RAW}/geo/GSE213982_series_matrix.txt.gz" || \
  echo "WARNING: Check exact matrix filename under https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213982/matrix/"

  # 2. Bryois 2022 Zenodo Record Metadata (lists all 8 cell-type files + md5)
  fetch "https://zenodo.org/api/records/5543734" "${OUT_RAW}/bryois_zenodo_5543734/zenodo_5543734_meta.json"
  python3 -c "
import json, pathlib
meta = json.loads(pathlib.Path('${OUT_RAW}/bryois_zenodo_5543734/zenodo_5543734_meta.json').read_text())
for f in meta.get('files', []):
    print(f['key'], round(f['size']/1024**2, 1), 'MB', f['links']['self'])
"
fi

echo "Done fetching automated manifests. For BrainMeta BESD & MetaBrain portal links, see 03_Aim2_Aim3_Server_Scripts/README.md."
