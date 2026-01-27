# ============================================
# STEP13J — Merge CytoTRACE2 scores into strict h5ad
# - reads strict cytotrace2-ready h5ad
# - reads scores CSV
# - writes scored h5ad to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

IN_H5AD_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/INPUTS"
SCORES_BASE  = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RESULTS"
OUT_BASE     = f"{ROOT}/outputs/STEP13_CYTOTRACE2/MERGED"
REPORT_DIR   = f"{ROOT}/outputs/STEP13_CYTOTRACE2/REPORTS"
os.makedirs(OUT_BASE, exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
rows = []

for s in samples:
    in_h5ad = f"{IN_H5AD_BASE}/{s}/spatial_{s}__strict__cytotrace2_ready.h5ad"
    score_csv = f"{SCORES_BASE}/{s}/cytotrace2_scores_{s}_human.csv"

    adata = sc.read_h5ad(in_h5ad)
    scores = pd.read_csv(score_csv, index_col=0)  # index = spot_id

    # Ensure index names match spot ids
    scores.index = scores.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)

    # Align scores to spots
    adata.obs["cytotrace2_potency_score"] = scores.loc[adata.obs_names, "cytotrace2_potency_score"].values

    # Save
    out_dir = f"{OUT_BASE}/{s}"
    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/spatial_{s}__strict__cytotrace2_scored.h5ad"
    adata.write_h5ad(out_path)

    # QC
    n_missing = int(pd.isna(adata.obs["cytotrace2_potency_score"]).sum())
    rows.append({
        "sample": s,
        "n_spots": adata.n_obs,
        "n_missing_scores": n_missing,
        "out_path": out_path,
        "saved_ok": os.path.exists(out_path)
    })

    print(f"\n{s}")
    print(" - spots:", adata.n_obs)
    print(" - missing scores:", n_missing)
    print(" - saved OK?:", os.path.exists(out_path))

qc = pd.DataFrame(rows)
qc_path = f"{REPORT_DIR}/merge_qc_scores_into_h5ad.csv"
qc.to_csv(qc_path, index=False)

print("\nSAVED TO DRIVE ✅ (merge QC)")
print("Path:", qc_path)
print("Saved OK?:", os.path.exists(qc_path))
print(qc)
