# SLC12A5 shared causal locus — analysis code

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![R 4.1+](https://img.shields.io/badge/R-4.1%2B-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reproducible analysis code for:

> **SLC12A5 is a shared causal locus for depression and anxiety, with opposing inherited and post-mortem effects on expression**

This repository contains **code only**. Raw data, derived binaries, figures and the manuscript are not included; acquisition and placement instructions are in [`data/README.md`](data/README.md).

---

## What the pipeline does

| Stage | Question | Method |
| :--- | :--- | :--- |
| 1. Locus preparation | Which variants enter fine-mapping? | 1 Mb dual-trait windows, allele alignment, 1000G EUR LD |
| 2. Multi-causal colocalization | Is one variant shared, or two linked variants? | `susieR::susie_rss` (L = 10) → `coloc::coloc.susie`; `coloc.abf` and `coloc::sensitivity` as controls |
| 3. Regulatory layer | Which cell type carries it, and in which direction? | Bryois 2022 eight-cell-type *cis*-eQTL colocalization, risk-allele harmonization, bulk-brain SMR/HEIDI |
| 4. Pseudobulk construction | How to aggregate 4.7 GB counts on 16 GB RAM? | Two-pass streaming C engine, zero full-matrix allocation |
| 5. Statistical modelling | Does the patient transcriptome follow the inherited direction? | `limma-voom` + `duplicateCorrelation`, subtype × diagnosis interaction, `limma::camera`, sex-stratified replication |
| 6. Composition and communication | Are the signals confounded? | Logit proportion models across 71 donors; per-donor LIANA+ (69 donors, 30,235 edges) |

```
dual-trait GWAS  ──►  fine-mapping  ──►  coloc.susie  ──►  shared vs distinct variant
                                    │
brain cell-type cis-eQTL ──► coloc ─┴─► allele harmonization ──► direction of effect
                                    │
GSE213982 counts ──► C streaming pseudobulk ──► limma-voom / camera / interaction
                                    │
                              ──► per-donor LIANA+ / composition models
```

---

## Repository layout

```
scripts/    23 files, named after the archived pipeline so any figure or table can be traced back
notebooks/  3 Colab Pro notebooks for the memory-heavy steps
tests/      syntax and hygiene checks (no data required)
data/       where downloaded inputs go (git-ignored)
```

| Stage | Script | Role |
| :--- | :--- | :--- |
| 0. cohort | `01_parse_geo_series_matrix_metadata.py` | Parse GEO series matrix into donor/diagnosis/sex table |
| 1. locus prep | `00_prepare_local_aim2_aligned_loci.py` | 1 Mb dual-trait windows, allele alignment, 1000G EUR LD |
| | `01_extract_gwas_loci_and_finemap.py` | Locus extraction and fine-mapping inputs (Aim 2 entry) |
| 2. colocalization | `04_run_aim2_susie_coloc_ctwas_hyprcoloc.R` | `susie_rss` + `coloc.susie`, single-causal controls |
| | `03_run_coloc_abf_when_eqtl_available.R` | ABF colocalization once eQTL inputs exist |
| | `12a_run_p0_coloc_abf_sensitivity.R` | Prior grid p12 ∈ [1e-6, 1e-4] + native sensitivity plot |
| 3. regulatory layer | `12b_harmonize_bryois_alleles.py` | Risk-allele harmonization (233/233 Bryois SNPs) |
| | `02_integrate_smr_heidi_twas_celltype.py` | SMR/HEIDI/TWAS × cell-type triangulation |
| | `05_audit_16panel_brain_smr_and_oxphos.py` | 16-panel bulk-brain SMR audit |
| | `04_run_smr_heidi_when_besd_available.sh` | SMR/HEIDI execution guard |
| | `03_server_fix_anx_freq_and_run_smr_multi.sh` | ANX frequency fix + multi-trait SMR (server) |
| 4. pseudobulk | `09_run_local_stream_pseudobulk.py` | Driver for the C engines |
| | `09a_stream_gse213982_pseudobulk.c` | Two-pass streaming aggregator (360.9 M entries, 374 MB peak) |
| | `09b_extract_cell_totals.c` | Library size / QC counter |
| 5. models | `10_run_p0_limma_camera_interaction.R` | `voom` + `duplicateCorrelation`, subtype interaction, `camera`, sex-stratified |
| | `11_run_p0_per_donor_liana.py` | Per-donor LIANA+ (69 donors, 30,235 edges) + OLS |
| 6. integration | `01_local_celltype_magma_spatial_integration.py` | Aim 1 scDRS × MAGMA × spatial integration |
| | `08_aim1_aim2_aim3_mechanistic_integration.py` | Cross-scale master tables |
| | `05_server_aim3_multilayer_integrate_and_resource_plan.py` | Multilayer integration and resource budget |
| 7. figures | `13_plot_publication_figures_ggplot2_pheatmap.R` | Figures 1–5 (native R) |
| | `07_plot_aim2_publication_figures.py` | Aim 2 colocalization figure |
| 8. resources | `02_download_aim2_aim3_public_resources.sh` | Public resource fetch with disk guardrail |
| | `02_local_ram_benchmark_probe.py` | Peak-memory probe used to size the C design |

Script numbering mirrors the archived working directory and is kept on purpose: files are not
renumbered, so any path or table reported in the manuscript can be traced back to the script
that produced it.

## Quick start

```bash
git clone https://github.com/chushuta-cmyk/slc12a5-shared-causal-locus.git
cd SLC12A5-shared-causal-locus

make c-engine          # compiles both C engines with -O3
make check             # syntax + repository hygiene tests, no data needed
```

To run the heavy stages, place the inputs as described in [`data/README.md`](data/README.md) and point the scripts at them:

```bash
export GSE213982_DIR=/path/to/GSE213982        # MatrixMarket counts + cell/gene tables
export MAITRA_ESM3=/path/to/maitra2023_esm3.xlsx
export PB_CACHE_DIR=/path/to/pseudobulk_cache  # optional override for stage 5

python3 scripts/09_run_local_stream_pseudobulk.py --check-only
python3 scripts/09_run_local_stream_pseudobulk.py
Rscript scripts/10_run_p0_limma_camera_interaction.R
```

---

## Environment actually used

| Component | Version |
| :--- | :--- |
| Python | 3.13.7 |
| scanpy / anndata / liana | 1.11.4 / 0.12.2 / 1.7.3 |
| numpy / pandas | 2.3.3 / 2.3.3 |
| R | 4.5.1 |
| limma / edgeR / coloc / susieR | 3.66.0 / 4.8.2 / 5.2.3 / 0.14.2 |
| Compiler | Apple clang 21 (C99, `-O3`) |

`liana` is required by stages 4 and 6 and is **not** a transitive dependency of the core scientific stack; install it explicitly (`pip install liana`). Environments that alias `python3` to a sandboxed interpreter without `liana` will fail at import — pin the interpreter you validated.

The full run was executed on a 16 GB Apple-silicon laptop: 55 s and 374 MB peak RSS to aggregate 360,887,743 non-zero entries into 2,431 donor × subtype pseudobulks.

---

## Data availability

No data are redistributed here. Upstream sources:

- **GSE213982** — human dorsolateral prefrontal cortex snRNA-seq (GEO). Cell and gene tables plus the MatrixMarket counts are needed for stage 4.
- **Bryois et al. 2022** — eight-brain-cell-type *cis*-eQTL summary statistics (Zenodo archive cited in the paper).
- **1000 Genomes Phase 3 EUR** — LD reference for fine-mapping (n = 503).
- **GWAS summary statistics** — depression and anxiety, European ancestry; sources are cited in the manuscript.

Derived caches (`*.bin`), downloaded archives, figures and the manuscript are excluded by `.gitignore` and must be regenerated.

---

## Citation

Citation details will be added on submission. Until then, please cite this repository by URL and the underlying datasets listed above.

## License

MIT — see [LICENSE](LICENSE).
