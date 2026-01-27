# ============================================
# STEP12 — Build STRICT cancer-rich sets
# Strategy:
# 1) Candidate spots already defined as Cancer_Epithelial > 0.10 (from STEP11)
# 2) Strict = top 10% (Q90) within candidates
# 3) Fallback: if strict too small -> take Top-N from candidates
# Saves ALL outputs to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Inputs
RESULTS_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/RESULTS"                 # full stereoscope outputs
CAND_BASE   = f"{ROOT}/outputs/STEP11_CANCER_RICH_SPOTS"                       # candidate h5ad (thr>0.10)

# # Outputs
OUT_BASE    = f"{ROOT}/outputs/STEP12_STRICT_CANCER_RICH"
os.makedirs(OUT_BASE, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# # Column name we created earlier (stereo_ + safe celltype name)
CANCER_COL = "stereo_cancer_epithelial"

# # Strict rule params
Q_STRICT = 0.90           # top 10%
MIN_STRICT_SPOTS = 200    # guarantee enough spots for downstream CytoTRACE2
TOPN_FALLBACK = 200       # fallback size (same as MIN_STRICT_SPOTS)

summary_rows = []
dist_rows = []

for s in samples:
    # ---------------------------
    # A) Distribution on FULL sample (for report + hocaya göstermek için)
    # ---------------------------
    full_h5ad = f"{RESULTS_DIR}/{s}/spatial_{s}_stereoscope_result.h5ad"
    ad_full = sc.read_h5ad(full_h5ad)

    if CANCER_COL not in ad_full.obs.columns:
        raise ValueError(f"{CANCER_COL} missing in full stereoscope result for {s}")

    x_full = ad_full.obs[CANCER_COL].astype(float).values
    dist_rows.append({
        "sample": s,
        "n_spots_full": ad_full.n_obs,
        "min": float(np.min(x_full)),
        "p50": float(np.quantile(x_full, 0.50)),
        "p75": float(np.quantile(x_full, 0.75)),
        "p90": float(np.quantile(x_full, 0.90)),
        "p95": float(np.quantile(x_full, 0.95)),
        "p99": float(np.quantile(x_full, 0.99)),
        "max": float(np.max(x_full)),
    })

    # ---------------------------
    # B) Load CANDIDATE set (thr>0.10) that we already saved in STEP11
    # ---------------------------
    cand_h5ad = f"{CAND_BASE}/{s}/spatial_{s}__cancer_rich_0.10.h5ad"
    ad_cand = sc.read_h5ad(cand_h5ad)

    if CANCER_COL not in ad_cand.obs.columns:
        raise ValueError(f"{CANCER_COL} missing in candidate h5ad for {s}")

    x = ad_cand.obs[CANCER_COL].astype(float)

    n_cand = ad_cand.n_obs

    # ---------------------------
    # C) Strict = Q90 within candidates
    # ---------------------------
    q_thr = float(x.quantile(Q_STRICT))
    mask_strict = (x >= q_thr)
    n_strict = int(mask_strict.sum())
    method_used = f"Q{int(Q_STRICT*100)}"

    # ---------------------------
    # D) Fallback: if too few strict spots, take Top-N from candidates
    # ---------------------------
    if n_cand == 0:
        # no candidates: save empty result (should not happen, but safe)
        ad_strict = ad_cand.copy()
        ad_strict = ad_strict[:0, :].copy()
        method_used = "EMPTY_CANDIDATES"
        q_thr = np.nan
    else:
        if (n_strict < MIN_STRICT_SPOTS) and (n_cand >= MIN_STRICT_SPOTS):
            # pick TOPN_FALLBACK spots by Cancer Epithelial proportion
            top_idx = x.sort_values(ascending=False).head(TOPN_FALLBACK).index
            ad_strict = ad_cand[top_idx, :].copy()
            method_used = f"TOP_{TOPN_FALLBACK}_FALLBACK"
            n_strict = ad_strict.n_obs
        else:
            # keep quantile-based strict (or candidate smaller than MIN -> that's still okay)
            ad_strict = ad_cand[mask_strict, :].copy()

    # ---------------------------
    # E) Save strict outputs to Drive
    # ---------------------------
    out_dir = f"{OUT_BASE}/{s}"
    os.makedirs(out_dir, exist_ok=True)

    # spot list
    spots_csv = f"{out_dir}/strict_cancer_rich_spots.csv"
    pd.DataFrame({"spot_id": ad_strict.obs_names.tolist()}).to_csv(spots_csv, index=False)

    # strict h5ad
    out_h5ad = f"{out_dir}/spatial_{s}__strict_cancer_rich.h5ad"
    ad_strict.uns["strict_selection"] = {
        "candidate_rule": f"{CANCER_COL} > 0.10",
        "strict_rule": f"top {int((1-Q_STRICT)*100)}% within candidates (Q{int(Q_STRICT*100)})",
        "quantile_threshold_value": q_thr,
        "min_strict_spots": MIN_STRICT_SPOTS,
        "fallback": f"top {TOPN_FALLBACK} if strict < {MIN_STRICT_SPOTS}",
        "method_used": method_used
    }
    ad_strict.write_h5ad(out_h5ad)

    print(f"\n{s}")
    print(" - candidates (thr>0.10):", n_cand)
    print(" - strict spots:", ad_strict.n_obs)
    print(" - method used:", method_used)
    print(" - saved CSV OK?:", os.path.exists(spots_csv))
    print(" - saved H5AD OK?:", os.path.exists(out_h5ad))

    summary_rows.append({
        "sample": s,
        "n_candidates_thr_0.10": n_cand,
        "q_strict": Q_STRICT,
        "q_threshold_value_in_candidates": q_thr,
        "method_used": method_used,
        "n_strict_spots": ad_strict.n_obs,
        "spots_csv_path": spots_csv,
        "strict_h5ad_path": out_h5ad
    })

# ---------------------------
# F) Save summary + distribution reports to Drive
# ---------------------------
summary_df = pd.DataFrame(summary_rows)
dist_df = pd.DataFrame(dist_rows)

summary_path = f"{OUT_BASE}/strict_cancer_rich_summary.csv"
dist_path    = f"{OUT_BASE}/cancer_epithelial_distribution_fullsample.csv"

summary_df.to_csv(summary_path, index=False)
dist_df.to_csv(dist_path, index=False)

print("\nSAVED TO DRIVE ✅ (summary)")
print("Path:", summary_path)
print("Saved OK?:", os.path.exists(summary_path))

print("\nSAVED TO DRIVE ✅ (distribution report)")
print("Path:", dist_path)
print("Saved OK?:", os.path.exists(dist_path))

print("\nDONE ✅")
