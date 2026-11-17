#!/usr/bin/env python3
"""P0 Per-Donor LIANA+ Verification (Biological Replicates N <= 71, Batch-Adjusted OLS)."""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-donor LIANA+ differential communication test")
    p.add_argument("--pb-dir", type=Path, default=SCRIPT_DIR / "outputs" / "local_pb_cache")
    p.add_argument("--out-tsv", type=Path, default=SCRIPT_DIR / "outputs" / "tables" / "p0_liana_per_donor_mdd_vs_ctrl_stats.tsv")
    p.add_argument("--check-only", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.check_only:
        print(f"P0_LIANA_PREFLIGHT: pb_dir_exists={int(args.pb_dir.exists())}, out_tsv={args.out_tsv.name}")
        return 0
    import scipy.sparse as sp
    import anndata as ad
    import liana as li
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests

    cell_map = np.fromfile(args.pb_dir / "cell_map.bin", dtype=np.int32).reshape(-1, 3)
    qc_orig_idx = np.fromfile(args.pb_dir / "qc_pass_orig_idx.bin", dtype=np.int32)
    tot_cnt = np.fromfile(args.pb_dir / "qc_pass_tot_cnt.bin", dtype=np.int64)
    lr_mask = np.fromfile(args.pb_dir / "lr_mask.bin", dtype=np.uint8)
    donors_df = pd.read_csv(args.pb_dir / "donor_metadata_complete.tsv", sep="\t")
    subs_df = pd.read_csv(args.pb_dir / "subtypes.tsv", sep="\t")
    genes_df = pd.read_csv(args.pb_dir / "genes.tsv", sep="\t")

    lr_orig = np.where(lr_mask == 1)[0]
    gene_to_lr = np.full(len(genes_df), -1, dtype=np.int32)
    gene_to_lr[lr_orig] = np.arange(len(lr_orig), dtype=np.int32)

    tri = np.fromfile(args.pb_dir / "lr_triplets.bin", dtype=np.int32).reshape(-1, 3)
    norm_vals = np.log1p(tri[:, 2].astype(np.float32) * (10000.0 / tot_cnt[tri[:, 1]].astype(np.float32)))
    X_csr = sp.csr_matrix((norm_vals, (tri[:, 1], gene_to_lr[tri[:, 0]])), shape=(len(qc_orig_idx), len(lr_orig)), dtype=np.float32)

    obs_df = pd.DataFrame({
        "donor_id": donors_df.set_index("donor_idx").loc[cell_map[qc_orig_idx, 0], "donor_id"].to_numpy(),
        "author_subtype": subs_df.set_index("subtype_idx").loc[cell_map[qc_orig_idx, 1], "author_subtype"].to_numpy(),
    }, index=[f"c_{i}" for i in range(len(qc_orig_idx))])
    adata = ad.AnnData(X=X_csr, obs=obs_df, var=pd.DataFrame(index=genes_df.iloc[lr_orig]["gene"].to_numpy()))

    focus = ["InN1_PV", "InN4_VIP", "InN5_SST", "ExN2_L23", "ExN3_L46", "ExN10_L46"]
    adata = adata[adata.obs["author_subtype"].isin(focus)].copy()

    per_donor_rows = []
    for did in sorted(donors_df["donor_id"].unique()):
        sub = adata[adata.obs["donor_id"] == did].copy()
        vc = sub.obs["author_subtype"].value_counts()
        keep_t = vc[vc >= 10].index.tolist()
        if not any(t.startswith("InN") for t in keep_t) or not any(t.startswith("ExN") for t in keep_t):
            continue
        sub = sub[sub.obs["author_subtype"].isin(keep_t)].copy()
        sub.obs["author_subtype"] = sub.obs["author_subtype"].astype(str).astype("category")
        li.mt.cellphonedb(sub, groupby="author_subtype", expr_prop=0.05, n_perms=10, use_raw=False, verbose=False, return_all_lrs=True)
        rd = sub.uns["liana_res"]
        rd = rd[rd["source"].astype(str).str.startswith("InN") & rd["target"].astype(str).str.startswith("ExN")].copy()
        rd["donor_id"] = did
        per_donor_rows.append(rd[["donor_id", "source", "target", "ligand_complex", "receptor_complex", "lr_means"]])

    df_all = pd.concat(per_donor_rows, ignore_index=True).merge(donors_df[["donor_id", "diagnosis", "Sex", "Batch", "Chemistry"]], on="donor_id", how="left")
    stats_rows = []
    for (src, tgt, lig, rec), grp in df_all.groupby(["source", "target", "ligand_complex", "receptor_complex"], observed=True):
        nc, nm = (grp["diagnosis"] == "Control").sum(), (grp["diagnosis"] == "MDD").sum()
        if nc < 8 or nm < 8 or (grp["lr_means"] > 0).sum() < 15:
            continue
        y = np.log1p(grp["lr_means"].to_numpy())
        d_bin = (grp["diagnosis"] == "MDD").astype(float).to_numpy()
        X = pd.concat([pd.Series(1.0, index=grp.index, name="const"), pd.Series(d_bin, index=grp.index, name="MDD"),
                       pd.get_dummies(grp["Batch"], drop_first=True, dtype=float)], axis=1)
        X = X.loc[:, (X.std(axis=0) > 0) | (X.columns == "const")]
        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            X = pd.concat([pd.Series(1.0, index=grp.index, name="const"), pd.Series(d_bin, index=grp.index, name="MDD"),
                           pd.get_dummies(grp["Sex"], drop_first=True, dtype=float), pd.get_dummies(grp["Chemistry"], drop_first=True, dtype=float)], axis=1)
            X = X.loc[:, (X.std(axis=0) > 0) | (X.columns == "const")]
        mod = sm.OLS(y, X).fit()
        ci = mod.conf_int().loc["MDD"]
        mc = float(grp.loc[grp["diagnosis"] == "Control", "lr_means"].mean())
        mm = float(grp.loc[grp["diagnosis"] == "MDD", "lr_means"].mean())
        stats_rows.append({
            "source": src, "target": tgt, "ligand": lig, "receptor": rec,
            "n_ctrl": int(nc), "n_mdd": int(nm), "mean_ctrl": mc, "mean_mdd": mm,
            "log2_fc": float(np.log2((mm + 1e-6) / (mc + 1e-6))),
            "ols_log1p_beta": float(mod.params["MDD"]), "se": float(mod.bse["MDD"]),
            "ci_low": float(ci[0]), "ci_high": float(ci[1]),
            "t_stat": float(mod.tvalues["MDD"]), "pvalue": float(mod.pvalues["MDD"]),
        })
    st = pd.DataFrame(stats_rows)
    st["FDR"] = multipletests(st["pvalue"], method="fdr_bh")[1]
    st.sort_values("pvalue").to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[SUCCESS] Per-donor LIANA complete: {(st['FDR'] < 0.05).sum()} / {len(st)} pairs FDR<0.05.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
