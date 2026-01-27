# ============================================
# STEP15B — Use pathology "classification" to define invasive FRONT vs CORE
# - reads pathology_annotation_candidates.csv from Drive
# - picks best file per sample automatically (max overlap, then label_col_score)
# - reads tissue_positions_list.csv to get grid coordinates
# - builds neighbor graph (k=6) using array_row/array_col
# - defines:
#   invasive_front = invasive spot with ANY neighbor in {stroma/adipo, normal ductal, DCIS}
#   tumor_core     = invasive spot with ALL neighbors invasive (within k)
# - merges with CytoTRACE2 scores (strict sets only)
# - saves per-sample spot table + stats + overall summary to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from scipy.stats import mannwhitneyu

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# Inputs
CAND_REPORT = f"{ROOT}/outputs/STEP15_INVASIVE_FRONT/REPORTS/pathology_annotation_candidates.csv"
MERGED_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/MERGED"  # strict scored h5ad per sample (has cytotrace score)
POS_BASE    = f"{ROOT}/extracted/spatial/spatial"        # tissue_positions_list.csv per sample

# Outputs
OUT_BASE = f"{ROOT}/outputs/STEP15_INVASIVE_FRONT"
os.makedirs(OUT_BASE, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
SCORE_COL = "cytotrace2_potency_score"

# ------------- helpers -------------
def read_positions(pos_path):
    # 10x legacy format (no header): barcode,in_tissue,array_row,array_col,pxl_row,pxl_col
    df = pd.read_csv(pos_path, header=None)
    df = df.iloc[:, :6]
    df.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    df["barcode"] = df["barcode"].astype(str)
    return df

def read_table(fp):
    # csv first, else tsv
    try:
        return pd.read_csv(fp)
    except Exception:
        return pd.read_csv(fp, sep="\t")

def normalize_class(x):
    return str(x).strip().lower()

def classify_pathology(val):
    """Map raw classification -> coarse category using keyword rules."""
    v = normalize_class(val)

    # invasive
    if "invasive" in v:
        return "invasive"

    # dcis
    if "dcis" in v:
        return "dcis"

    # stroma / adipose
    if ("stroma" in v) or ("adipo" in v):
        return "stroma_adipose"

    # normal ductal/epithelial
    if ("normal" in v) or ("duct" in v):
        return "normal_ductal"

    # lymphocyte aggregates
    if ("lymph" in v) or ("aggregate" in v):
        return "lymph_agg"

    return "other"

# ------------- load candidate report -------------
rep = pd.read_csv(CAND_REPORT)

# Keep only CANDIDATE rows
rep = rep[rep["status"] == "CANDIDATE"].copy()

# Choose ONE best file per sample:
# sort by overlap_spots desc, then label_col_score desc
rep = rep.sort_values(["sample", "overlap_spots", "label_col_score"], ascending=[True, False, False])

best = rep.groupby("sample").head(1).copy()

print("Best pathology file per sample:")
print(best[["sample", "file_path", "barcode_col", "suggested_label_col", "overlap_spots", "label_col_score"]])

all_stats = []
all_rows = []

for s in samples:
    # paths
    h5ad_path = f"{MERGED_BASE}/{s}/spatial_{s}__strict__cytotrace2_scored.h5ad"
    pos_path  = f"{POS_BASE}/{s}_spatial/tissue_positions_list.csv"

    if not os.path.exists(h5ad_path):
        print(f"{s}: missing scored h5ad, skip.")
        continue
    if not os.path.exists(pos_path):
        print(f"{s}: missing tissue_positions_list.csv, skip.")
        continue

    # select pathology file info
    if s not in best["sample"].values:
        print(f"{s}: no pathology candidate found in report, skip.")
        continue

    row = best[best["sample"] == s].iloc[0]
    meta_fp = row["file_path"]
    bc_col = row["barcode_col"]
    label_col = row["suggested_label_col"]  # expected: classification

    # 1) load positions (all spots grid)
    pos = read_positions(pos_path)
    pos = pos.set_index("barcode")

    # 2) load pathology table
    meta = read_table(meta_fp)
    meta[bc_col] = meta[bc_col].astype(str)
    meta = meta.set_index(bc_col)

    if label_col not in meta.columns:
        print(f"{s}: label_col '{label_col}' not in {meta_fp}, skip.")
        continue

    # 3) merge positions + pathology label
    df = pos.join(meta[[label_col]], how="left")
    df.rename(columns={label_col: "pathology_label_raw"}, inplace=True)

    # 4) map coarse categories
    df["pathology_category"] = df["pathology_label_raw"].apply(classify_pathology)

    # 5) build neighbor graph on array coordinates (k=6)
    coords = df[["array_col", "array_row"]].astype(float).values
    nn = NearestNeighbors(n_neighbors=7, metric="euclidean")  # 7 includes self
    nn.fit(coords)
    neigh_idx = nn.kneighbors(return_distance=False)

    # helper to get neighbor categories
    cats = df["pathology_category"].values
    barcodes = df.index.values

    # 6) define front/core on ALL invasive spots
    is_invasive = (cats == "invasive")
    front = np.zeros(len(df), dtype=bool)
    core  = np.zeros(len(df), dtype=bool)

    target_contact = set(["stroma_adipose", "normal_ductal", "dcis"])

    for i in range(len(df)):
        if not is_invasive[i]:
            continue
        neigh = neigh_idx[i, 1:]  # exclude self
        neigh_cats = set(cats[neigh])
        # front: touches any non-invasive target
        if len(neigh_cats.intersection(target_contact)) > 0:
            front[i] = True
        # core: all neighbors invasive (no contact with anything else)
        if all(cats[j] == "invasive" for j in neigh):
            core[i] = True

    df["is_invasive"] = is_invasive
    df["is_front"] = front
    df["is_core"] = core
    df["region"] = np.where(front, "front", np.where(core, "core", "invasive_other"))

    # 7) merge CytoTRACE2 scores (strict spots only)
    ad = sc.read_h5ad(h5ad_path)
    ad.obs_names = ad.obs_names.astype(str)
    scores = ad.obs[[SCORE_COL]].copy()
    scores.index.name = "barcode"

    df = df.join(scores, how="left")  # scores only for strict spots

    # keep only spots that have CytoTRACE2 score
    df_scored = df[~df[SCORE_COL].isna()].copy()

    # stats: compare front vs core using Mann-Whitney U
    front_scores = df_scored.loc[df_scored["region"] == "front", SCORE_COL].astype(float).values
    core_scores  = df_scored.loc[df_scored["region"] == "core",  SCORE_COL].astype(float).values

    stat = {
        "sample": s,
        "pathology_file": meta_fp,
        "positions_file": pos_path,
        "n_scored_spots": int(df_scored.shape[0]),
        "n_front_scored": int(len(front_scores)),
        "n_core_scored": int(len(core_scores)),
        "front_median": float(np.median(front_scores)) if len(front_scores) else np.nan,
        "core_median": float(np.median(core_scores)) if len(core_scores) else np.nan,
        "front_mean": float(np.mean(front_scores)) if len(front_scores) else np.nan,
        "core_mean": float(np.mean(core_scores)) if len(core_scores) else np.nan,
        "mw_u_pvalue": np.nan
    }

    if len(front_scores) >= 10 and len(core_scores) >= 10:
        # two-sided test
        _, p = mannwhitneyu(front_scores, core_scores, alternative="two-sided")
        stat["mw_u_pvalue"] = float(p)

    all_stats.append(stat)

    # save per-sample spot table
    out_table = f"{OUT_BASE}/front_core_spot_table_{s}.csv"
    df_scored[[
        "in_tissue", "array_row", "array_col",
        "pathology_label_raw", "pathology_category",
        "is_invasive", "is_front", "is_core", "region",
        SCORE_COL
    ]].to_csv(out_table, index=True)

    print(f"\n{s}: saved spot table ✅", os.path.exists(out_table))
    print("  scored spots:", stat["n_scored_spots"],
          "| front:", stat["n_front_scored"],
          "| core:", stat["n_core_scored"],
          "| p:", stat["mw_u_pvalue"])

# save stats
stats_df = pd.DataFrame(all_stats)
stats_path = f"{OUT_BASE}/per_sample_front_core_stats.csv"
stats_df.to_csv(stats_path, index=False)

# overall summary: simple aggregation (median of medians)
summary = {
    "n_samples": int(stats_df.shape[0]),
    "median_front_median": float(np.nanmedian(stats_df["front_median"])) if len(stats_df) else np.nan,
    "median_core_median": float(np.nanmedian(stats_df["core_median"])) if len(stats_df) else np.nan
}
summary_df = pd.DataFrame([summary])
summary_path = f"{OUT_BASE}/front_core_summary_all_samples.csv"
summary_df.to_csv(summary_path, index=False)

print("\nSAVED TO DRIVE ✅ (per-sample stats):", stats_path, os.path.exists(stats_path))
print("SAVED TO DRIVE ✅ (overall summary):", summary_path, os.path.exists(summary_path))
print("\nStats preview:")
print(stats_df)
