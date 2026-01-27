# ============================================
# STEP13F — Build CPM (no-log) input + run CytoTRACE2 (HUMAN) for 1142243F
# Saves input + results to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import pickle

from cytotrace2_py.cytotrace2_py import cytotrace2

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Start from cytotrace2-ready h5ad (has layers['counts'])
IN_H5AD = f"{ROOT}/outputs/STEP13_CYTOTRACE2/INPUTS/1142243F/spatial_1142243F__strict__cytotrace2_ready.h5ad"

OUT_IN_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RUN_INPUTS/1142243F"
OUT_RES_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RESULTS/1142243F"
os.makedirs(OUT_IN_DIR, exist_ok=True)
os.makedirs(OUT_RES_DIR, exist_ok=True)

IN_TXT = f"{OUT_IN_DIR}/cytotrace2_input_1142243F_CPM.txt"
OUT_CSV = f"{OUT_RES_DIR}/cytotrace2_scores_1142243F_human.csv"
OUT_PKL = f"{OUT_RES_DIR}/cytotrace2_output_full_1142243F_human.pkl"

# ------------------------
# 1) Build CPM (no log) matrix from counts
# ------------------------
adata = sc.read_h5ad(IN_H5AD)

X = adata.layers["counts"]
if issparse(X):
    X = X.toarray()
X = X.astype(np.float64)

# CPM normalize per spot (row)
row_sums = X.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1.0
X_cpm = (X / row_sums) * 1e6  # CPM (no log)

expr = pd.DataFrame(X_cpm.T, index=adata.var_names, columns=adata.obs_names)  # genes x spots
expr.to_csv(IN_TXT, sep="\t")

print("SAVED TO DRIVE ✅ (CPM input)")
print("Path:", IN_TXT)
print("Saved OK?:", os.path.exists(IN_TXT))
print("Matrix shape (genes x spots):", expr.shape)

# ------------------------
# 2) Run CytoTRACE2 with correct species = HUMAN
# ------------------------
print("\nRunning CytoTRACE2 (species='human') ...")
res = cytotrace2(input_path=IN_TXT, species="human", plot=False, verbose=True)

print("\nCytoTRACE2 finished ✅")

# Save full output
with open(OUT_PKL, "wb") as f:
    pickle.dump(res, f)
print("Saved full output ✅:", OUT_PKL, "OK?", os.path.exists(OUT_PKL))

# Extract scores (handle dict/attr)
score_series = None
if isinstance(res, dict):
    for k in ["potency_score", "cytotrace2_score", "score", "scores"]:
        if k in res:
            score_series = res[k]
            break

if score_series is None:
    for attr in ["potency_score", "cytotrace2_score", "score", "scores"]:
        if hasattr(res, attr):
            score_series = getattr(res, attr)
            break

if score_series is None:
    if isinstance(res, dict):
        print("Keys:", list(res.keys()))
    raise ValueError("Could not extract CytoTRACE2 scores from output.")

scores = pd.Series(score_series)
scores_df = scores.to_frame(name="cytotrace2_potency_score")
scores_df.index.name = "spot_id"
scores_df.to_csv(OUT_CSV)

print("\nSAVED TO DRIVE ✅ (scores CSV)")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("Preview:")
print(scores_df.head())
