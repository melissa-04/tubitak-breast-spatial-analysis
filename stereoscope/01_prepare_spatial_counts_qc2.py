# ============================================
# STEP09 — Prepare SPATIAL QC2 for Stereoscope
# - load normalized files
# - set X = raw counts (from layers["counts"])
# - convert counts to int
# - save stereoscope-ready h5ad to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import scanpy as sc

PROJECT_ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Input: normalized QC2 outputs (these have layers["counts"] + layers["log2_cpm"])
IN_DIR  = f"{PROJECT_ROOT}/outputs/STEP08_NORMALIZATION/QC2"

# # Output: stereoscope-ready QC2 (X = integer counts)
OUT_DIR = f"{PROJECT_ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SPATIAL_QC2_COUNTS"
os.makedirs(OUT_DIR, exist_ok=True)

print("Input dir :", IN_DIR)
print("Output dir:", OUT_DIR)

# # Find files we created in STEP08
in_files = sorted([f for f in os.listdir(IN_DIR) if f.endswith("__log2cpm.h5ad")])
print("\nFound normalized QC2 files:", len(in_files))
for f in in_files:
    print(" -", f)

for f in in_files:
    in_path  = f"{IN_DIR}/{f}"
    out_path = f"{OUT_DIR}/{f.replace('__log2cpm.h5ad', '__counts_for_stereoscope.h5ad')}"

    adata = sc.read_h5ad(in_path)

    # # Safety checks
    if "counts" not in adata.layers:
        raise ValueError(f"'counts' layer missing in: {in_path}")

    # # Set X = raw counts (Stereoscope input)
    X = adata.layers["counts"]

    # # counts float32 geldiği için int'e çevir (round güvenli çünkü orijinalde integer idi)
    # # (sparse ise de çalışır)
    if hasattr(X, "A"):  # sparse matrix
        X_int = X.copy()
        X_int.data = np.rint(X_int.data).astype(np.int64)
        adata.X = X_int
    else:  # dense
        adata.X = np.rint(X).astype(np.int64)

    # # Küçük özet
    total_counts = int(adata.X.sum())
    print(f"\nSaved prep for: {f}")
    print(" - shape:", adata.shape)
    print(" - total counts (X.sum):", total_counts)

    # # Save
    adata.write_h5ad(out_path)
    print(" - saved OK?:", os.path.exists(out_path))
    print(" - out:", out_path)

print("\nDONE: Spatial QC2 stereoscope-prep completed.")
