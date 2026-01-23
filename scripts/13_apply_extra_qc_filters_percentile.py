import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

in_root  = base / "outputs" / "STEP07_H5AD_MAIN"
out_root = base / "outputs" / "STEP07B_H5AD_QC_FILTERED"
out_root.mkdir(parents=True, exist_ok=True)

# percentile settings
LOW_Q  = 0.01
HIGH_Q = 0.99
MT_MAX = 20.0  # percent

rows = []

for s in samples:
    in_h5ad = in_root / s / f"{s}__MAIN_keep_in_tissue_no_artefact.h5ad"
    adata = sc.read_h5ad(in_h5ad)

    # ensure qc metrics exist; compute mt% if possible
    mt_mask = adata.var_names.str.upper().str.startswith("MT-")
    if mt_mask.any():
        adata.var["mt"] = mt_mask
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
        has_mt = True
    else:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        has_mt = False

    # thresholds (per sample)
    tc = adata.obs["total_counts"].to_numpy()
    ng = adata.obs["n_genes_by_counts"].to_numpy()

    tc_lo, tc_hi = np.quantile(tc, [LOW_Q, HIGH_Q])
    ng_lo, ng_hi = np.quantile(ng, [LOW_Q, HIGH_Q])

    mask = (tc >= tc_lo) & (tc <= tc_hi) & (ng >= ng_lo) & (ng <= ng_hi)

    if has_mt:
        mask = mask & (adata.obs["pct_counts_mt"].to_numpy() <= MT_MAX)

    before = adata.n_obs
    adata_f = adata[mask, :].copy()
    after = adata_f.n_obs
    dropped = before - after

    # save
    out_dir = out_root / s
    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = out_dir / f"{s}__QC2_filtered.h5ad"
    adata_f.write_h5ad(out_h5ad)

    row = {
        "sample": s,
        "before_spots": int(before),
        "after_spots": int(after),
        "dropped_spots": int(dropped),
        "dropped_%": float(100 * dropped / max(1, before)),
        "tc_low_q01": float(tc_lo),
        "tc_high_q99": float(tc_hi),
        "ng_low_q01": float(ng_lo),
        "ng_high_q99": float(ng_hi),
        "mt_filter_used": bool(has_mt),
        "mt_max_%": float(MT_MAX if has_mt else np.nan),
    }
    rows.append(row)

    print(f"✅ {s}: before={before} -> after={after} (dropped={dropped}, {row['dropped_%']:.2f}%) | saved: {out_h5ad.name}")

# combined report
df = pd.DataFrame(rows).sort_values("sample")
report_path = out_root / "STEP07B__QC2_filter_report_ALL.csv"
df.to_csv(report_path, index=False)

print("\n✅ Saved report:", report_path)
display(df)
