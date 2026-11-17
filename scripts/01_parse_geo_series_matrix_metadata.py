#!/usr/bin/env python3
"""
Parse Official NCBI GEO Series Matrix (`*_series_matrix.txt.gz`) into Clean Metadata TSV
======================================================================================
Extracts `!Sample_geo_accession`, `!Sample_title`, and all `!Sample_characteristics_ch1`
key:value fields into a normalized donor metadata TSV ready for `03_colab_pro_upstream_snRNA_spatial_scDRS.ipynb`
or any GEO bulk/snRNA validation cohort.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path
import pandas as pd


def normalize_col_name(raw: str) -> str:
    s = re.sub(r"[^0-9a-zA-Z]+", "_", raw.strip().lower()).strip("_")
    return s or "attr"


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def parse_geo_series_matrix(matrix_path: Path) -> pd.DataFrame:
    sample_titles: list[str] = []
    sample_acc: list[str] = []
    char_rows: list[list[str]] = []

    with open_text(matrix_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            tag = row[0]
            if tag == "!series_matrix_table_begin":
                break
            elif tag == "!Sample_title":
                sample_titles = [c.strip() for c in row[1:]]
            elif tag == "!Sample_geo_accession":
                sample_acc = [c.strip() for c in row[1:]]
            elif tag == "!Sample_characteristics_ch1":
                char_rows.append([c.strip() for c in row[1:]])

    if not sample_acc:
        raise ValueError(f"No !Sample_geo_accession header found in {matrix_path}")

    n_samples = len(sample_acc)
    records: list[dict[str, str]] = []
    for i in range(n_samples):
        rec: dict[str, str] = {
            "geo_accession": sample_acc[i],
            "sample_title": sample_titles[i] if i < len(sample_titles) else "",
        }
        # Infer short donor/sample prefix if formatted like "F1" or "F1_dlPFC"
        m_id = re.match(r"^([FM]\d+(?:_\d+)?)\b", rec["sample_title"])
        if m_id:
            sid = m_id.group(1)
            rec["sample_id"] = sid
            rec["donor_id"] = "M24" if sid == "M24_2" else sid

        for r_idx, ch in enumerate(char_rows):
            if i >= len(ch) or not ch[i]:
                continue
            val = ch[i]
            if ":" in val:
                k, v = val.split(":", 1)
                cname = normalize_col_name(k)
                if cname in rec and rec[cname] != v.strip():
                    cname = f"{cname}_{r_idx}"
                rec[cname] = v.strip()
            else:
                rec[f"characteristic_{r_idx}"] = val

        # Standardize diagnosis column if common phenotype labels exist
        for pheno_col in ("diagnosis", "disease_state", "disorder", "condition", "phenotype", "best_overall_response"):
            if pheno_col in rec and "diagnosis" not in rec:
                rec["diagnosis"] = rec[pheno_col]
        records.append(rec)

    return pd.DataFrame(records)


def main(argv: list[str] | None = None) -> pd.DataFrame:
    p = argparse.ArgumentParser(description="Parse NCBI GEO *_series_matrix.txt.gz clinical metadata.")
    p.add_argument("--matrix", type=Path, required=True, help="Path to *_series_matrix.txt.gz")
    p.add_argument("--out", type=Path, required=True, help="Output TSV path")
    args = p.parse_args(argv)

    df = parse_geo_series_matrix(args.matrix)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Parsed {len(df)} samples with {len(df.columns)} metadata fields -> {args.out}")
    print("Columns:", df.columns.tolist())
    return df


if __name__ == "__main__":
    main()
