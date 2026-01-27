# ============================================
# STEP15A — Find pathology annotation files + best label columns (per sample)
# - scans extracted/ for csv/tsv/txt
# - matches files to each sample by barcode overlap
# - suggests which column looks like pathology label
# Saves report to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import numpy as np
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# strict scored h5ad (spot IDs here)
SCORED_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/MERGED"

# search only extracted (fast + relevant)
SEARCH_ROOTS = [
    f"{ROOT}/extracted/metadata",
    f"{ROOT}/extracted/spatial",
    f"{ROOT}/extracted"  # fallback
]

OUT_DIR = f"{ROOT}/outputs/STEP15_INVASIVE_FRONT/REPORTS"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_CSV = f"{OUT_DIR}/pathology_annotation_candidates.csv"

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# keywords to recognize pathology-like columns/values
col_kw = ["path", "annot", "region", "hist", "morph", "class", "label", "compart", "tissue"]
val_kw = ["invasive", "dcis", "stroma", "adipo", "normal", "duct", "lymph", "aggregate"]

def try_read_table(path):
    # try csv then tsv
    try:
        df = pd.read_csv(path)
        return df
    except Exception:
        try:
            df = pd.read_csv(path, sep="\t")
            return df
        except Exception:
            return None

def find_barcode_col(df):
    # common barcode column names
    for c in df.columns:
        cl = str(c).lower()
        if cl in ["barcode", "barcodes", "spot_id", "spotid", "id", "unnamed: 0"]:
            return c
    # heuristic: a column with many strings ending with "-1"
    for c in df.columns:
        if df[c].dtype == object:
            vals = df[c].astype(str).head(200).values
            if np.mean([v.endswith("-1") for v in vals]) > 0.6:
                return c
    return None

def score_label_column(series):
    # score by whether values look like pathology categories
    s = series.dropna().astype(str).str.lower()
    if len(s) == 0:
        return 0
    hits = 0
    for k in val_kw:
        if s.str.contains(k).any():
            hits += 1
    return hits

# collect candidate data files
candidate_files = []
for sr in SEARCH_ROOTS:
    if not os.path.exists(sr):
        continue
    for r, d, files in os.walk(sr):
        for f in files:
            lf = f.lower()
            if lf.endswith((".csv", ".tsv", ".txt")):
                # keep likely metadata-ish names
                if any(k in lf for k in ["meta", "annot", "path", "region", "hist", "label", "morph", "class"]):
                    candidate_files.append(os.path.join(r, f))

candidate_files = sorted(list(set(candidate_files)))
print("Candidate files found:", len(candidate_files))

rows = []

for s in samples:
    scored_path = f"{SCORED_BASE}/{s}/spatial_{s}__strict__cytotrace2_scored.h5ad"
    ad = sc.read_h5ad(scored_path)
    spot_ids = set(ad.obs_names.astype(str).tolist())

    best_hits = []

    for fp in candidate_files:
        df = try_read_table(fp)
        if df is None or df.shape[0] == 0:
            continue

        bc_col = find_barcode_col(df)
        if bc_col is None:
            continue

        barcodes = set(df[bc_col].astype(str).tolist())
        overlap = len(spot_ids.intersection(barcodes))
        if overlap == 0:
            continue

        # candidate label columns by name
        label_cols = [c for c in df.columns if any(k in str(c).lower() for k in col_kw)]
        # score each label col by value keywords
        best_col = ""
        best_col_score = -1
        best_col_nuniq = None

        for c in label_cols:
            scv = score_label_column(df[c])
            if scv > best_col_score:
                best_col_score = scv
                best_col = c
                best_col_nuniq = int(df[c].nunique(dropna=True))

        best_hits.append((overlap, fp, bc_col, best_col, best_col_score, best_col_nuniq, df.shape[0], df.shape[1]))

    # keep top 3 by overlap
    best_hits = sorted(best_hits, key=lambda x: x[0], reverse=True)[:3]

    if len(best_hits) == 0:
        rows.append({
            "sample": s,
            "status": "NO_MATCH_FOUND",
            "file_path": "",
            "barcode_col": "",
            "overlap_spots": 0,
            "suggested_label_col": "",
            "label_col_score": 0,
            "label_n_unique": "",
            "n_rows": "",
            "n_cols": ""
        })
        continue

    for (overlap, fp, bc_col, best_col, best_col_score, best_col_nuniq, nrows, ncols) in best_hits:
        rows.append({
            "sample": s,
            "status": "CANDIDATE",
            "file_path": fp,
            "barcode_col": bc_col,
            "overlap_spots": overlap,
            "suggested_label_col": best_col,
            "label_col_score": best_col_score,
            "label_n_unique": best_col_nuniq,
            "n_rows": nrows,
            "n_cols": ncols
        })

report = pd.DataFrame(rows)
report.to_csv(OUT_CSV, index=False)

print("\nSAVED TO DRIVE ✅")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("\nPreview:")
print(report.head(20))
