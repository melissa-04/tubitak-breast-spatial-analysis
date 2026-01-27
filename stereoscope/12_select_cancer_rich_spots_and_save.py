# ============================================
# STEP11 — Select cancer-rich spots using Stereoscope proportions
# - thresholds: >0.10 and >0.50 for Cancer Epithelial
# - saves spot lists (CSV) + filtered h5ad to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

RESULTS_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/RESULTS"
OUT_BASE    = f"{ROOT}/outputs/STEP11_CANCER_RICH_SPOTS"
os.makedirs(OUT_BASE, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# # Stereoscope'ta bu kolon adını ürettik: stereo_cancer_epithelial
CANCER_COL = "stereo_cancer_epithelial"

thresholds = [0.10, 0.50]
summary = []

for s in samples:
    in_h5ad = f"{RESULTS_DIR}/{s}/spatial_{s}_stereoscope_result.h5ad"
    adata = sc.read_h5ad(in_h5ad)

    if CANCER_COL not in adata.obs.columns:
        raise ValueError(f"{CANCER_COL} missing in {s}. Check stereoscope outputs.")

    out_dir = f"{OUT_BASE}/{s}"
    os.makedirs(out_dir, exist_ok=True)

    for thr in thresholds:
        # # Spot seçimi
        mask = adata.obs[CANCER_COL] > thr
        ad_sub = adata[mask, :].copy()

        # # Spot listesi (barcode/index)
        spot_ids = ad_sub.obs_names.to_list()
        spot_csv = f"{out_dir}/cancer_rich_spots_threshold_{thr:.2f}.csv"
        pd.DataFrame({"spot_id": spot_ids}).to_csv(spot_csv, index=False)

        # # Filtrelenmiş h5ad kaydı
        out_h5ad = f"{out_dir}/spatial_{s}__cancer_rich_{thr:.2f}.h5ad"
        ad_sub.write_h5ad(out_h5ad)

        # # Özet satırı
        summary.append({
            "sample": s,
            "threshold": thr,
            "n_spots_selected": ad_sub.n_obs,
            "spots_csv_path": spot_csv,
            "h5ad_path": out_h5ad
        })

        print(f"\n{s} | thr>{thr:.2f}")
        print("Selected spots:", ad_sub.n_obs)
        print("Saved CSV OK?:", os.path.exists(spot_csv))
        print("Saved H5AD OK?:", os.path.exists(out_h5ad))

# # Özet dosyasını Drive'a kaydet
summary_df = pd.DataFrame(summary)
summary_path = f"{OUT_BASE}/summary_cancer_rich_counts.csv"
summary_df.to_csv(summary_path, index=False)

print("\nSAVED TO DRIVE ✅ (summary)")
print("Path:", summary_path)
print("Saved OK?:", os.path.exists(summary_path))
print(summary_df)
