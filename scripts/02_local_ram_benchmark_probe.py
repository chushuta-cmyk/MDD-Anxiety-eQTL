#!/usr/bin/env python3
"""Safe Local Hardware & RAM Scaling Probe for Scheme 1 Upstream Workflows.

Runs each bounded AnnData size (5,000 / 10,000 / 20,000 cells x 20,000 genes) in an
isolated fresh Python subprocess so pymalloc arenas from smaller runs never mask
true incremental RSS growth. Measures both process RSS and active Scanpy memory
footprint, then scales to full GSE213982 (~160,000 nuclei x ~36,000 genes + Scrublet
doublets + dense scaled HVG + Harmony + 6 GB snRNA reference Cell2location).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent.parent
BM_DIR = BASE_DIR / "outputs" / "benchmarks"
BM_DIR.mkdir(parents=True, exist_ok=True)


def run_single_worker(n_cells: int, n_genes: int = 20000, density: float = 0.06) -> dict:
    """Execute one bounded sparse AnnData preprocessing pipeline in the current clean process."""
    import anndata as ad
    import psutil
    import scanpy as sc
    import scipy.sparse as sp
    from sklearn.neighbors import NearestNeighbors

    proc = psutil.Process(os.getpid())
    rss_init = proc.memory_info().rss / (1024 ** 2)
    t0 = time.perf_counter()

    rng = np.random.default_rng(42)
    nnz = int(n_cells * n_genes * density)
    rows = rng.integers(0, n_cells, size=nnz, dtype=np.int32)
    cols = rng.integers(0, n_genes, size=nnz, dtype=np.int32)
    vals = rng.poisson(2.0, size=nnz).astype(np.float32) + 1.0
    X = sp.csr_matrix((vals, (rows, cols)), shape=(n_cells, n_genes), dtype=np.float32)
    del rows, cols, vals

    adata = ad.AnnData(X=X)
    adata.obs["batch"] = [f"B{i % 20}" for i in range(n_cells)]
    rss_loaded = proc.memory_info().rss / (1024 ** 2)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, svd_solver="arpack", n_comps=50)
    NearestNeighbors(n_neighbors=15, n_jobs=1).fit(adata_hvg.obsm["X_pca"]).kneighbors(adata_hvg.obsm["X_pca"])
    rss_peak = proc.memory_info().rss / (1024 ** 2)
    elapsed = time.perf_counter() - t0

    return {
        "n_cells": n_cells,
        "n_genes": n_genes,
        "rss_init_mb": round(rss_init, 2),
        "rss_loaded_mb": round(rss_loaded, 2),
        "rss_post_pca_knn_mb": round(rss_peak, 2),
        "net_workload_rss_mb": round(rss_peak - rss_init, 2),
        "elapsed_sec": round(elapsed, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-cells", type=int, default=0)
    args = parser.parse_args()

    if args.worker_cells > 0:
        res = run_single_worker(args.worker_cells)
        print(json.dumps(res))
        return

    total_sys_ram_gb = int(subprocess.check_output(["sysctl", "-n", "hw.memsize"]).decode().strip()) / (1024 ** 3)
    steps = [5000, 10000, 20000]
    results = []
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"

    for n in steps:
        p = subprocess.run(
            ["/usr/bin/time", "-l", sys.executable, __file__, "--worker-cells", str(n)],
            capture_output=True, text=True, env=env, check=True,
        )
        # Parse JSON from last line of stdout
        row = json.loads(p.stdout.strip().splitlines()[-1])
        # Extract macOS peak memory footprint & max RSS from /usr/bin/time -l stderr
        max_rss_bytes = None
        peak_footprint_bytes = None
        for line in p.stderr.splitlines():
            if "maximum resident set size" in line:
                max_rss_bytes = int(line.strip().split()[0])
            elif "peak memory footprint" in line:
                peak_footprint_bytes = int(line.strip().split()[0])
        row["macos_max_rss_mb"] = round(max_rss_bytes / (1024 ** 2), 2) if max_rss_bytes else row["rss_post_pca_knn_mb"]
        row["macos_peak_footprint_mb"] = round(peak_footprint_bytes / (1024 ** 2), 2) if peak_footprint_bytes else None
        results.append(row)

    df = pd.DataFrame(results)
    # Use macOS peak_footprint_mb minus 5k-offset baseline for honest physical memory demand
    y_col = "macos_peak_footprint_mb" if df["macos_peak_footprint_mb"].notna().all() else "macos_max_rss_mb"
    slope_mb_per_cell, _ = np.polyfit(df["n_cells"], df[y_col], 1)

    # Full GSE213982 upstream (~160k nuclei x 36k genes -> 1.8x gene width)
    # Plus Scrublet simulated doublets + PyTorch/Cell2location 6GB reference + Harmony
    proj_pca_only_gb = (slope_mb_per_cell * 160000 * 1.8) / 1024.0
    proj_full_upstream_gb = proj_pca_only_gb * 1.65  # conservative multiplier for Scrublet + Harmony + Cell2location

    summary = {
        "local_machine_ram_gb": round(total_sys_ram_gb, 2),
        "safe_local_headroom_gb": round(total_sys_ram_gb * 0.70, 2),
        "measured_isolated_subprocesses": results,
        "measured_slope_mb_per_1k_cells_20k_genes": round(float(slope_mb_per_cell * 1000.0), 2),
        "projected_160k_36k_qc_pca_gb": round(float(proj_pca_only_gb), 2),
        "projected_full_upstream_with_scrublet_harmony_c2l_gb": round(float(proj_full_upstream_gb), 2),
        "local_feasible_upstream": bool(proj_full_upstream_gb < (total_sys_ram_gb * 0.70)),
        "recommended_runtime": (
            "Google Colab Pro (High-RAM 51-54 GB + A100/L4 GPU) for upstream GSE213982/Cell2location/scDRS; "
            "Local Mac (16 GB) for downstream pseudobulk/MAGMA/spatial integration (< 0.4 GB RAM)."
        ),
    }

    df.to_csv(BM_DIR / "local_ram_scaling_curve.tsv", sep="\t", index=False)
    with open(BM_DIR / "local_ram_probe.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
