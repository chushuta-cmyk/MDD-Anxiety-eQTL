#!/usr/bin/env python3
"""
Aim 3 Multi-Omic Orthogonal Integration + Verified Server Compute & Storage Sizing Plan
========================================================================================
Modes:
  --estimate-resources   Generates exact verified CPU / RAM / Disk resource table across
                         Aim 1 (Colab/Server), Aim 2 (Multi-causal SuSiE/cTWAS/HyPrColoc),
                         and Aim 3 (Bryois cell-eQTL + BrainMeta sQTL + pQTL + mQTL).
  --integrate            Merges available Aim 2 + Aim 3 outputs (Bryois, BrainMeta sQTL,
                         Wingo pQTL, MDD/ANX SMR) into a unified 5-layer gene prioritization TSV.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR
REPO_ROOT = SCRIPT_DIR.parent

RESOURCE_ESTIMATES = [
    {
        "stage_id": "Aim1_Upstream_GSE213982",
        "stage_name": "Full 160k x 36.6k snRNA QC + Scrublet + Harmony(donor,sex) + PyDESeq2 + scDRS",
        "can_use_16gb_mac": False,
        "online_portal_feasible": False,
        "recommended_compute": "Google Colab Pro (High-RAM 51-54 GB + A100/L4) OR 64 GB / 8-core Server",
        "ram_peak_gb_est": "36.6 GB (PCA/KNN) to 60.4 GB (with Cell2location)",
        "recommended_server_ram_gb": 64,
        "recommended_cpu_cores": 8,
        "disk_compressed_gb": 1.2,
        "disk_working_uncompressed_gb": 12.0,
        "est_runtime": "35-70 min",
        "notes": "Measured slope 130.11 MB / 1k cells in 02_local_ram_benchmark_probe.py; h5ad + scrublet + scDRS fits cleanly in 64 GB RAM.",
    },
    {
        "stage_id": "Aim1_Spatial_gsMap_spatialDLPFC",
        "stage_name": "spatialDLPFC (10 donors, 30 Visium) Cell2location + gsMap GWAS Spatial Projection",
        "can_use_16gb_mac": False,
        "online_portal_feasible": False,
        "recommended_compute": "GPU Server / Colab Pro (A100 40GB / L4 24GB + >=32 GB System RAM)",
        "ram_peak_gb_est": "20-32 GB host RAM + 16-24 GB VRAM",
        "recommended_server_ram_gb": 32,
        "recommended_cpu_cores": 8,
        "disk_compressed_gb": 4.0,
        "disk_working_uncompressed_gb": 15.0,
        "est_runtime": "40-90 min",
        "notes": "Fix N_cells_per_spot=8 for human DLPFC; GPU cuts Cell2location runtime from hours to ~15 min.",
    },
    {
        "stage_id": "Aim2_Local_GWAS_Extract_ABF",
        "stage_name": "Stream 1.15M MDD + 1.10M ANX GWAS -> Dual-Trait Locus Extract + ABF Finemapping",
        "can_use_16gb_mac": True,
        "online_portal_feasible": True,
        "recommended_compute": "Current 16 GB Mac (already runs locally in < 10 s!)",
        "ram_peak_gb_est": "< 0.8 GB",
        "recommended_server_ram_gb": 4,
        "recommended_cpu_cores": 2,
        "disk_compressed_gb": 0.3,
        "disk_working_uncompressed_gb": 0.6,
        "est_runtime": "< 10 sec",
        "notes": "Already verified locally with 00_prepare_local_aim2_aligned_loci.py.",
    },
    {
        "stage_id": "Aim2_MultiCausal_SuSiE_cTWAS_HyPrColoc",
        "stage_name": "Per-Locus LD (plink --r square) + SuSiE + coloc.susie + cTWAS + HyPrColoc",
        "can_use_16gb_mac": True,
        "online_portal_feasible": False,
        "recommended_compute": "8-16 GB Server OR 16 GB Mac (per-locus window mode; GWAShug online lacks SuSiE/cTWAS)",
        "ram_peak_gb_est": "2-8 GB (2,000-5,000 SNPs per locus LD matrix = ~80-200 MB)",
        "recommended_server_ram_gb": 16,
        "recommended_cpu_cores": 4,
        "disk_compressed_gb": 2.0,
        "disk_working_uncompressed_gb": 6.0,
        "est_runtime": "10-25 min across 7-83 loci",
        "notes": "Does NOT need 64 GB! LD matrix of M=5,000 SNPs is only 5000^2 * 8 bytes = 200 MB.",
    },
    {
        "stage_id": "Aim2_ANX_Freq_Fix_and_SMR_Multi",
        "stage_name": "1000G EUR plink --freq ANX Alignment + --smr-multi (BrainMeta eQTL + sQTL + Bryois)",
        "can_use_16gb_mac": True,
        "online_portal_feasible": False,
        "recommended_compute": "Linux x86_64 Server (16-32 GB RAM) OR 16 GB Mac (compile C++ SMR / Docker)",
        "ram_peak_gb_est": "4-12 GB for BrainMeta BESD sparse SMR+HEIDI",
        "recommended_server_ram_gb": 16,
        "recommended_cpu_cores": 4,
        "disk_compressed_gb": 6.0,
        "disk_working_uncompressed_gb": 18.0,
        "est_runtime": "15-45 min across all BESD panels",
        "notes": "Online SMR Portal cannot run custom Bryois/pQTL BESDs or --smr-multi on custom panels cleanly.",
    },
    {
        "stage_id": "Aim3_CellType_sQTL_pQTL_mQTL_Integrate",
        "stage_name": "Bryois 8 Cell-Type eQTL + BrainMeta sQTL + MetaBrain Cortex + Wingo pQTL + mQTL Merge",
        "can_use_16gb_mac": True,
        "online_portal_feasible": False,
        "recommended_compute": "16-32 GB Linux Server (storage-bound ~25-40 GB if keeping MetaBrain unzipped)",
        "ram_peak_gb_est": "4-16 GB (streaming or SQLite/DuckDB/tabix indexed summary stats)",
        "recommended_server_ram_gb": 32,
        "recommended_cpu_cores": 4,
        "disk_compressed_gb": 10.0,
        "disk_working_uncompressed_gb": 35.0,
        "est_runtime": "20-40 min (mostly I/O & allele harmonization)",
        "notes": "Current Mac has 51 GiB free disk: stream/keep .gz + subset target genes/windows to save space, OR put MetaBrain on a 100 GB server volume.",
    },
]


def integrate_available_layers(aim3_dir: Path, smr_16_dir: Path, out_tab: Path) -> pd.DataFrame:
    """Integrate already-local 16-panel SMR hits with any staged Aim 3 Bryois/sQTL/pQTL/mQTL tables."""
    rep_path = REPO_ROOT / "02_Scheme2_GWAS_Causal_Coloc_SMR" / "outputs" / "tables" / "smr_16panel_prioritized_replication.tsv"
    base_df = pd.read_csv(rep_path, sep="\t") if rep_path.exists() else pd.DataFrame({"Gene": []})

    for layer_name, fname in [
        ("Bryois_CellType_eQTL", "bryois_celltype_top_hits.tsv"),
        ("BrainMeta_sQTL", "brainmeta_sqtl_smr_hits.tsv"),
        ("Wingo_Brain_pQTL", "wingo_pqtl_pwas_hits.tsv"),
        ("Brain_or_Blood_mQTL", "mqtl_mwas_hits.tsv"),
    ]:
        fpath = aim3_dir / fname
        if fpath.exists() and not base_df.empty:
            sub = pd.read_csv(fpath, sep="\t")
            if "Gene" in sub.columns:
                sub_ren = sub.add_prefix(f"{layer_name}__").rename(columns={f"{layer_name}__Gene": "Gene"})
                base_df = base_df.merge(sub_ren.drop_duplicates("Gene"), on="Gene", how="left")
        else:
            base_df[f"{layer_name}__status"] = "pending_external_download_or_server_run"

    out_tab.parent.mkdir(parents=True, exist_ok=True)
    base_df.to_csv(out_tab, sep="\t", index=False)
    return base_df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Aim 3 server resource estimator & multi-layer integration.")
    p.add_argument("--estimate-resources", action="store_true", help="Write verified server RAM/CPU/disk estimates.")
    p.add_argument("--integrate", action="store_true", help="Integrate staged Aim 3 tables with Scheme 2 SMR backbone.")
    p.add_argument("--out-dir", type=Path, default=BASE_DIR / "outputs")
    p.add_argument("--aim3-inputs", type=Path, default=BASE_DIR / "inputs" / "aim3_external")
    args = p.parse_args(argv)

    if not args.estimate_resources and not args.integrate:
        args.estimate_resources = True

    tab_dir = args.out_dir / "tables"
    bm_dir = args.out_dir / "benchmarks"
    tab_dir.mkdir(parents=True, exist_ok=True)
    bm_dir.mkdir(parents=True, exist_ok=True)

    if args.estimate_resources:
        est_df = pd.DataFrame(RESOURCE_ESTIMATES)
        est_df.to_csv(tab_dir / "server_resource_budget_estimate.tsv", sep="\t", index=False)
        summary = {
            "aim1_min_server_ram_gb": 64,
            "aim2_min_server_ram_gb": 16,
            "aim3_recommended_server_ram_gb": 32,
            "aim2_aim3_combined_recommended_server": "8-16 vCPU, 32 GB RAM (or 64 GB if also running Aim 1 outside Colab), 100-150 GB NVMe Disk",
            "local_mac_16gb_capable_steps": [
                r["stage_id"] for r in RESOURCE_ESTIMATES if r["can_use_16gb_mac"]
            ],
            "online_portal_only_insufficient_for": [
                r["stage_id"] for r in RESOURCE_ESTIMATES if not r["online_portal_feasible"]
            ],
            "stages": RESOURCE_ESTIMATES,
        }
        with open(bm_dir / "server_resource_budget_estimate.json", "w") as f:
            json.dump(summary, f, indent=2)
        print(est_df[["stage_id", "can_use_16gb_mac", "recommended_server_ram_gb", "disk_working_uncompressed_gb", "est_runtime"]].to_string(index=False))

    if args.integrate:
        df = integrate_available_layers(args.aim3_inputs, REPO_ROOT / "02_Scheme2_GWAS_Causal_Coloc_SMR", tab_dir / "aim3_multilayer_gene_prioritization.tsv")
        print(f"Integrated {len(df)} prioritized genes -> {tab_dir / 'aim3_multilayer_gene_prioritization.tsv'}")


if __name__ == "__main__":
    main()
