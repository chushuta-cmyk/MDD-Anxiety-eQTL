# Where the inputs go

Nothing in `data/` is tracked by git. Create the following layout, or override the
paths with the environment variables listed below.

```
data/
├── GSE213982/
│   ├── GSE213982_combined_counts_matrix.mtx                  (~4.7 GB)
│   ├── GSE213982_combined_counts_matrix_cells_columns.csv.gz
│   └── GSE213982_combined_counts_matrix_genes_rows.csv.gz
├── maitra2023_esm3.xlsx        # donor / batch metadata, SupplementaryData1
└── pb_cache/                   # written by the C engines, read by stage 5
```

| Variable | Used by | Default |
| :--- | :--- | :--- |
| `GSE213982_DIR` | `scripts/09_run_local_stream_pseudobulk.py` | `data/GSE213982` |
| `MAITRA_ESM3` | same script | `data/maitra2023_esm3.xlsx` |
| `PB_CACHE_DIR` | `scripts/10_run_p0_limma_camera_interaction.R` | `data/pb_cache` |

Sources:

- **GSE213982** — GEO series with the processed MatrixMarket counts. The `.mtx` file must be in
  MatrixMarket coordinate form; the streaming engines read it twice and never hold it in memory.
- **Maitra et al. 2023, SupplementaryData1** — donor-level diagnosis, sex, batch, chemistry,
  sequencing platform and QC metrics. Twelve batch levels, which is what the models adjust for.
- **Bryois et al. 2022 eQTL summary statistics** — not stored here; stage 3 and the notebooks
  read them from the Zenodo archive cited in the paper.
- **1000 Genomes Phase 3 EUR** — LD reference for fine-mapping; stage 1 expects per-locus
  genotype panels.

Expected scale if you want to confirm that aggregation ran correctly: 360,887,743 non-zero
entries, 121,284 QC-passing nuclei, 2,431 non-empty donor × subtype pseudobulks, 339 MB of
pseudobulk counts.
