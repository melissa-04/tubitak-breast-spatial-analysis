# ============================================
# STEP14-FINAL — Plot using sample-specific tissue_positions file
# - chooses positions file that contains the sample id in its path
# - regenerates PNG with "_fixed" suffix
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

ROOT = "/content/drive/MyDrive/TÜBİTAK"
KEY = "cytotrace2_potency_score"

SCORED_WITHOUT_COORDS = f"{ROOT}/outputs/STEP13_CYTOTRACE2/MERGED"
OUT_DIR = f"{ROOT}/outputs/STEP14_FIGURES/CYTOTRACE2_SPATIAL"
os.makedirs(OUT_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# find all tissue_positions files
all_pos = []
for r, d, files in os.walk(f"{ROOT}/extracted"):
    for f in files:
        if "tissue_positions" in f.lower() and f.lower().endswith(".csv"):
            all_pos.append(os.path.join(r, f))
all_pos = sorted(list(set(all_pos)))

def read_positions(path):
    try:
        df = pd.read_csv(path)
        if "barcode" in df.columns:
            return df
    except Exception:
        pass
    df = pd.read_csv(path, header=None)
    df = df.iloc[:, :6]
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    return df

for s in samples:
    h5ad = f"{SCORED_WITHOUT_COORDS}/{s}/spatial_{s}__strict__cytotrace2_scored.h5ad"
    ad = sc.read_h5ad(h5ad)
    ad.obs_names = ad.obs_names.astype(str)

    # choose positions file that contains sample id in its path
    candidates = [p for p in all_pos if s in p]
    if len(candidates) == 0:
        print(f"{s}: no sample-specific tissue_positions found, skip.")
        continue

    pos_path = candidates[0]  # usually unique
    dfp = read_positions(pos_path)
    dfp["barcode"] = dfp["barcode"].astype(str)
    dfp = dfp.set_index("barcode")

    # choose coords columns
    if "pxl_row_in_fullres" in dfp.columns and "pxl_col_in_fullres" in dfp.columns:
        ycol, xcol = "pxl_row_in_fullres", "pxl_col_in_fullres"
    else:
        ycol, xcol = "pxl_row", "pxl_col"

    # align
    keep = [bc for bc in ad.obs_names if bc in dfp.index]
    ad = ad[keep, :].copy()
    coords = np.vstack([dfp.loc[bc, [xcol, ycol]].astype(float).values for bc in ad.obs_names])

    # plot
    out_png = f"{OUT_DIR}/cytotrace2_scatter_{s}_fixed.png"
    plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1], c=ad.obs[KEY].values, s=16)
    plt.gca().invert_yaxis()
    plt.title(f"{s} — CytoTRACE2 potency (strict) [fixed coords]")
    plt.colorbar()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close()

    print(f"{s}: saved FIXED PNG OK?:", os.path.exists(out_png), "| positions:", pos_path)

print("\nDONE ✅")
