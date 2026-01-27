# ============================================
# STEP14 — Attach coords from 10x tissue_positions*.csv + plot CytoTRACE2
# - finds best-matching tissue_positions file by barcode overlap
# - adds obsm['spatial'] to scored h5ad
# - saves PNG scatter + new h5ad to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

ROOT = "/content/drive/MyDrive/TÜBİTAK"

SCORED_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/MERGED"
OUT_DIR = f"{ROOT}/outputs/STEP14_FIGURES/CYTOTRACE2_SPATIAL"
OUT_H5AD_DIR = f"{OUT_DIR}/WITH_COORDS"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(OUT_H5AD_DIR, exist_ok=True)

KEY = "cytotrace2_potency_score"
samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# 1) Find all tissue_positions files under TÜBİTAK (anywhere)
pos_files = []
for r, d, files in os.walk(ROOT):
    for f in files:
        lf = f.lower()
        if lf in ["tissue_positions_list.csv", "tissue_positions.csv"] or "tissue_positions" in lf:
            pos_files.append(os.path.join(r, f))

pos_files = sorted(list(set(pos_files)))
print("tissue_positions files found:", len(pos_files))
for p in pos_files[:15]:
    print(" -", p)
if len(pos_files) > 15:
    print("(showing first 15)")

def read_positions(path):
    """Read tissue_positions file (header/no-header variants)."""
    # Try header-based first
    try:
        df = pd.read_csv(path)
        if "barcode" in df.columns:
            return df
    except Exception:
        pass

    # No-header legacy format
    df = pd.read_csv(path, header=None)
    # legacy columns: barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
    if df.shape[1] >= 6:
        df = df.iloc[:, :6]
        df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    return df

def best_positions_file(barcodes_set):
    """Pick the tissue_positions file with max barcode overlap."""
    best_path, best_overlap = None, -1
    for p in pos_files:
        try:
            df = read_positions(p)
            if "barcode" not in df.columns:
                continue
            bc = set(df["barcode"].astype(str).values)
            ov = len(barcodes_set.intersection(bc))
            if ov > best_overlap:
                best_overlap = ov
                best_path = p
        except Exception:
            continue
    return best_path, best_overlap

report_rows = []

for s in samples:
    scored_path = f"{SCORED_BASE}/{s}/spatial_{s}__strict__cytotrace2_scored.h5ad"
    if not os.path.exists(scored_path):
        print(f"{s}: scored h5ad not found, skip.")
        continue

    ad = sc.read_h5ad(scored_path)
    ad.obs_names = ad.obs_names.astype(str)
    barcodes = set(ad.obs_names.values)

    # 2) Find best tissue_positions file for this sample
    best_path, ov = best_positions_file(barcodes)

    if best_path is None or ov == 0:
        print(f"{s}: NO matching tissue_positions file found (overlap=0).")
        report_rows.append({"sample": s, "positions_file": "", "overlap_barcodes": ov, "status": "NO_MATCH"})
        continue

    dfp = read_positions(best_path)
    dfp["barcode"] = dfp["barcode"].astype(str)
    dfp = dfp.set_index("barcode")

    # 3) Choose coordinate columns
    # New format: pxl_row_in_fullres / pxl_col_in_fullres
    if "pxl_row_in_fullres" in dfp.columns and "pxl_col_in_fullres" in dfp.columns:
        ycol, xcol = "pxl_row_in_fullres", "pxl_col_in_fullres"
    # Legacy format: pxl_row / pxl_col
    elif "pxl_row" in dfp.columns and "pxl_col" in dfp.columns:
        ycol, xcol = "pxl_row", "pxl_col"
    # Fallback: array grid coordinates
    elif "array_row" in dfp.columns and "array_col" in dfp.columns:
        ycol, xcol = "array_row", "array_col"
    else:
        print(f"{s}: positions file has no usable coord columns.")
        report_rows.append({"sample": s, "positions_file": best_path, "overlap_barcodes": ov, "status": "NO_COORD_COLS"})
        continue

    # 4) Build coords aligned to strict spots
    keep = [bc for bc in ad.obs_names if bc in dfp.index]
    missing = len(ad.obs_names) - len(keep)
    if missing > 0:
        ad = ad[keep, :].copy()

    coords = np.vstack([dfp.loc[bc, [xcol, ycol]].astype(float).values for bc in ad.obs_names])
    ad.obsm["spatial"] = coords  # (x, y)

    # 5) Save new h5ad with coords
    out_h5ad = f"{OUT_H5AD_DIR}/spatial_{s}__strict__cytotrace2_scored_with_coords.h5ad"
    ad.write_h5ad(out_h5ad)

    # 6) Plot scatter (robust, no sc.pl.spatial)
    out_png = f"{OUT_DIR}/cytotrace2_scatter_{s}.png"
    plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1], c=ad.obs[KEY].values, s=16)
    plt.gca().invert_yaxis()
    plt.title(f"{s} — CytoTRACE2 potency (strict)")
    plt.colorbar()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close()

    print(f"{s}: saved PNG OK?: {os.path.exists(out_png)} | saved H5AD OK?: {os.path.exists(out_h5ad)}")
    print(f"   positions file used: {best_path} | overlap={ov} | dropped_missing={missing}")

    report_rows.append({
        "sample": s,
        "positions_file": best_path,
        "overlap_barcodes": ov,
        "dropped_missing_barcodes_in_strict": missing,
        "coord_columns_used": f"{xcol},{ycol}",
        "png_path": out_png,
        "h5ad_with_coords_path": out_h5ad,
        "status": "OK"
    })

# 7) Save report
rep = pd.DataFrame(report_rows)
rep_path = f"{OUT_DIR}/coords_source_report.csv"
rep.to_csv(rep_path, index=False)

print("\nSAVED TO DRIVE ✅ (coords source report)")
print("Path:", rep_path)
print("Saved OK?:", os.path.exists(rep_path))
print("\nDONE ✅")
