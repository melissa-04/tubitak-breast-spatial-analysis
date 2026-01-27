# ============================================
# STEP10F — QC check for Stereoscope outputs
# - verifies files exist
# - loads each h5ad
# - checks stereo_* columns and sum-to-one
# Saves a QC summary CSV to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"
RESULTS_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/RESULTS"
REPORT_DIR  = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/REPORTS"
os.makedirs(REPORT_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

rows = []
for s in samples:
    h5ad_path = f"{RESULTS_DIR}/{s}/spatial_{s}_stereoscope_result.h5ad"
    csv_path  = f"{RESULTS_DIR}/{s}/stereoscope_proportions_{s}.csv"

    exist_h5ad = os.path.exists(h5ad_path)
    exist_csv  = os.path.exists(csv_path)

    info = {
        "sample": s,
        "h5ad_exists": exist_h5ad,
        "csv_exists": exist_csv,
        "n_spots": None,
        "n_celltypes": None,
        "sum_min": None,
        "sum_max": None,
        "sum_mean": None,
        "n_nan_sums": None
    }

    if exist_h5ad:
        adata = sc.read_h5ad(h5ad_path)

        # stereo_* kolonlarını bul
        stereo_cols = [c for c in adata.obs.columns if str(c).startswith("stereo_")]

        if len(stereo_cols) > 0:
            mat = adata.obs[stereo_cols].to_numpy(dtype=float)
            sums = np.sum(mat, axis=1)

            info["n_spots"] = adata.n_obs
            info["n_celltypes"] = len(stereo_cols)
            info["sum_min"] = float(np.nanmin(sums))
            info["sum_max"] = float(np.nanmax(sums))
            info["sum_mean"] = float(np.nanmean(sums))
            info["n_nan_sums"] = int(np.isnan(sums).sum())

    rows.append(info)

qc = pd.DataFrame(rows)
out_path = f"{REPORT_DIR}/stereoscope_output_qc_summary.csv"
qc.to_csv(out_path, index=False)

print("SAVED TO DRIVE ✅")
print("Path:", out_path)
print("Saved OK?:", os.path.exists(out_path))
print("\nQC Summary:")
print(qc)
