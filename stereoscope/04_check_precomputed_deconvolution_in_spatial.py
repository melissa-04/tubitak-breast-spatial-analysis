# ============================================
# STEP09F — Check if spatial files already contain deconvolution columns
# Saves a CSV report to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import numpy as np
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Spatial stereoscope-ready counts dosyaların
SPATIAL_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SPATIAL_QC2_COUNTS"

# # scRNA'dan bildiğimiz major cell type isimleri (referans için)
SCRNA_REF = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE/GSE176078_scRNA_reference_counts__label_celltype_major.h5ad"
scrna = sc.read_h5ad(SCRNA_REF)
known_types = sorted(scrna.obs["cell_type"].unique().tolist())

# # Raporun kaydedileceği yer (Drive)
REPORT_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/REPORTS"
os.makedirs(REPORT_DIR, exist_ok=True)
OUT_CSV = f"{REPORT_DIR}/precomputed_deconvolution_check.csv"

print("Known cell types (from scRNA major):", known_types)
print("Spatial dir:", SPATIAL_DIR)

# # Spatial dosyaları listele
files = sorted([f for f in os.listdir(SPATIAL_DIR) if f.endswith("__counts_for_stereoscope.h5ad")])
print("Found spatial files:", len(files))

# # Kolon adlarında arayacağımız anahtar kelimeler
keywords = [
    "stereo", "deconv", "proportion", "fraction",
    "epithelial", "cancer", "t-cell", "tcell", "myeloid",
    "endothelial", "caf", "plasma", "b-cell", "bcell", "pvl"
]

rows = []

for f in files:
    path = f"{SPATIAL_DIR}/{f}"

    # # sadece obs'a bakacağız; backed okuma daha hafif
    adata = sc.read_h5ad(path, backed="r")
    obs = adata.obs

    # # aday kolonlar: ismi keyword içeriyor mu?
    name_hits = []
    for c in obs.columns:
        cl = str(c).lower()
        if any(k in cl for k in keywords):
            name_hits.append(c)

    # # aday kolonlar: sayısal mı ve 0-1 aralığına benziyor mu?
    numeric_hits = []
    for c in obs.columns:
        if pd.api.types.is_numeric_dtype(obs[c]):
            col = obs[c].dropna()
            if len(col) > 0:
                if (col.min() >= -0.05) and (col.max() <= 1.05):
                    numeric_hits.append(c)

    # # ikisinin kesişimi en güçlü aday
    strong = sorted(set(name_hits).intersection(set(numeric_hits)))

    rows.append({
        "spatial_file": f,
        "n_obs_columns": obs.shape[1],
        "name_hit_cols": ";".join(map(str, name_hits[:50])),
        "strong_candidate_cols_0to1": ";".join(map(str, strong[:50])),
        "n_strong_candidates": len(strong)
    })

    print(f"\n{f}")
    print(" - obs columns:", obs.shape[1])
    print(" - strong candidates (0-1 + keyword):", len(strong))

# # Raporu Drive'a kaydet
df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)

print("\nSAVED TO DRIVE ✅")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("\nPreview:")
print(df[["spatial_file", "n_strong_candidates"]])
