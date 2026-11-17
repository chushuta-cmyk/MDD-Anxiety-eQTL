#!/usr/bin/env python3
"""Compile & run low-memory C streaming pseudobulk aggregator on raw GSE213982 MatrixMarket counts.

Produces:
  - Donor x Subtype pseudobulk count tensor (pb_counts.bin + pb_groups_meta.tsv)
  - Sparse LIANA+ consensus LR single-cell triplets (lr_triplets.bin + qc_pass_tot_cnt.bin)
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[0]
# Point $GSE213982_DIR at the unpacked GEO supplementary files, or place them in ./data/GSE213982
DEFAULT_RAW_DIR = Path(os.environ.get("GSE213982_DIR", str(ROOT / "data" / "GSE213982")))
DEFAULT_ESM3 = Path(os.environ.get("MAITRA_ESM3", str(ROOT / "data" / "maitra2023_esm3.xlsx")))
DEFAULT_OUT = SCRIPT_DIR / "outputs" / "local_pb_cache"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Low-RAM C streaming GSE213982 pseudobulk aggregator")
    p.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    p.add_argument("--esm3-xlsx", type=Path, default=DEFAULT_ESM3)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    p.add_argument("--check-only", action="store_true", help="Compile C helpers and check paths only.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    c_pb = SCRIPT_DIR / "09a_stream_gse213982_pseudobulk.c"
    c_tot = SCRIPT_DIR / "09b_extract_cell_totals.c"
    assert c_pb.exists() and c_tot.exists()

    if args.check_only:
        mtx_ok = int((args.raw_dir / "GSE213982_combined_counts_matrix.mtx").exists())
        cc_ok = subprocess.run(["which", "cc"], capture_output=True).returncode == 0
        print(f"LOCAL_STREAM_PREFLIGHT: cc={int(cc_ok)}, raw_mtx={mtx_ok}, c_sources=2")
        return 0

    import liana as li
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cells_gz = args.raw_dir / "GSE213982_combined_counts_matrix_cells_columns.csv.gz"
    genes_gz = args.raw_dir / "GSE213982_combined_counts_matrix_genes_rows.csv.gz"
    mtx_file = args.raw_dir / "GSE213982_combined_counts_matrix.mtx"

    cells = pd.read_csv(cells_gz).iloc[:, 0].astype(str).str.strip('"')
    genes = pd.read_csv(genes_gz).iloc[:, 0].astype(str).str.strip('"')
    toks = cells.str.split(".", expand=True)
    donors_raw = toks[0].replace({"M24_2": "M24"})
    subtypes_raw = toks[3]

    uniq_donors = sorted(donors_raw.unique())
    uniq_subs = sorted(subtypes_raw.unique())
    d_map = {d: i for i, d in enumerate(uniq_donors)}
    s_map = {s: i for i, s in enumerate(uniq_subs)}

    cell_map = np.zeros((len(cells), 3), dtype=np.int32)
    cell_map[:, 0] = donors_raw.map(d_map).to_numpy()
    cell_map[:, 1] = subtypes_raw.map(s_map).to_numpy()
    cell_map[:, 2] = subtypes_raw.str.contains("Mix", case=False).astype(np.int32).to_numpy()
    cell_map.tofile(args.out_dir / "cell_map.bin")

    don_df = pd.DataFrame({"donor_idx": range(len(uniq_donors)), "donor_id": uniq_donors})
    sub_df = pd.DataFrame({"subtype_idx": range(len(uniq_subs)), "author_subtype": uniq_subs})
    don_df.to_csv(args.out_dir / "donors.tsv", sep="\t", index=False)
    sub_df.to_csv(args.out_dir / "subtypes.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_idx": range(len(genes)), "gene": genes}).to_csv(args.out_dir / "genes.tsv", sep="\t", index=False)

    res = li.resource.select_resource("consensus")
    lr_syms = set(res["ligand"].astype(str)).union(set(res["receptor"].astype(str)))
    lr_syms = lr_syms.union({x for g in lr_syms for x in g.split("_")})
    lr_mask = genes.isin(lr_syms).astype(np.uint8).to_numpy()
    lr_mask.tofile(args.out_dir / "lr_mask.bin")

    bin_pb = args.out_dir / "stream_pb"
    bin_tot = args.out_dir / "extract_totals"
    subprocess.run(["cc", "-O3", str(c_pb), "-o", str(bin_pb)], check=True)
    subprocess.run(["cc", "-O3", str(c_tot), "-o", str(bin_tot)], check=True)
    subprocess.run([
        str(bin_pb), str(mtx_file), str(args.out_dir / "cell_map.bin"),
        str(args.out_dir / "lr_mask.bin"), str(len(uniq_donors)), str(len(uniq_subs)), str(args.out_dir),
    ], check=True)
    subprocess.run([
        str(bin_tot), str(mtx_file), str(args.out_dir / "qc_pass_orig_idx.bin"),
        str(args.out_dir / "qc_pass_tot_cnt.bin"),
    ], check=True)

    # Merge Maitra 2023 SupplementaryData1 metadata
    if args.esm3_xlsx.exists():
        s1 = pd.read_excel(args.esm3_xlsx, sheet_name="SupplementaryData1", skiprows=2)
        s1 = s1[s1["Sample"].notna()].copy()
        s1["donor_id"] = s1["Sample"].astype(str).replace({"M24_2": "M24"})
        s1_u = s1.drop_duplicates("donor_id").copy()
        s1_u["diagnosis"] = s1_u["Condition"].map({"Case": "MDD", "Control": "Control"})
        don_full = don_df.merge(
            s1_u[["donor_id", "diagnosis", "Sex", "Batch", "Chemistry", "Sequencing",
                  "MedianGene_filtered", "MedianUMI_filtered", "SequencingSaturation", "NumCells_filtered"]],
            on="donor_id", how="left",
        )
        don_full.to_csv(args.out_dir / "donor_metadata_complete.tsv", sep="\t", index=False)
        grp = pd.read_csv(args.out_dir / "pb_groups.tsv", sep="\t")
        grp_full = grp.merge(don_full, on="donor_idx", how="left").merge(sub_df, on="subtype_idx", how="left")
        grp_full.to_csv(args.out_dir / "pb_groups_meta.tsv", sep="\t", index=False)

    print(f"[SUCCESS] Local streaming aggregation finished cleanly in {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
