# ============================================
# STEP16 — Run CytoTRACE2 on CANDIDATE sets
# for ALL samples + front vs core statistics
# ============================================

from google.colab import drive
drive.mount("/content/drive")

import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from cytotrace2_py import cytotrace2

ROOT = "/content/drive/MyDrive/TÜBİTAK"

CAND_BASE = f"{ROOT}/outputs/STEP15_INVASIVE_FRONT/CANDIDATE_H5AD"
LABEL_BASE = f"{ROOT}/outputs/STEP15_INVASIVE_FRONT/LABELS"
OUT_BASE = f"{ROOT}/outputs/STEP16_CANDIDATE_CYTOTRACE2"

os.makedirs(f"{OUT_BASE}/RESULTS", exist_ok=True)
os.makedirs(f"{OUT_BASE}/REPORTS", exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

all_rows = []

for s in samples:
    print(f"\n=== {s} ===")

    h5ad_path = f"{CAND_BASE}/{s}/candidate_{s}.h5ad"
    label_path = f"{LABEL_BASE}/{s}_front_core_labels.csv"

    ad = sc.read_h5ad(h5ad_path)

    # CPM
    X = ad.layers["counts"].astype(float)
    X_cpm = X / X.sum(axis=0) * 1e6
    ad.layers["cpm"] = X_cpm

    # write CytoTRACE2 input
    out_dir = f"{OUT_BASE}/RESULTS/{s}"
    os.makedirs(out_dir, exist_ok=True)
    in_txt = f"{out_dir}/cytotrace2_input_{s}_CPM.txt"

    pd.DataFrame(
        ad.layers["cpm"],
        index=ad.var_names,
        columns=ad.obs_names
    ).to_csv(in_txt, sep="\t")

    # run CytoTRACE2
    cytotrace2(
        input_path=in_txt,
        species="human",
        output_dir=out_dir,
        plot=False,
        verbose=False
    )

    # read scores
    score_csv = [f for f in os.listdir(out_dir) if f.endswith(".csv")][0]
    scores = pd.read_csv(f"{out_dir}/{score_csv}", index_col=0)

    # merge with front/core labels
    lab = pd.read_csv(label_path, index_col=0)
    df = lab.join(scores, how="inner")

    front = df[df["region"] == "front"]["cytotrace2_potency_score"]
    core  = df[df["region"] == "core"]["cytotrace2_potency_score"]

    if len(front) > 0 and len(core) > 0:
        u, p = mannwhitneyu(front, core, alternative="two-sided")
    else:
        p = np.nan

    row = {
        "sample": s,
        "n_front": len(front),
        "n_core": len(core),
        "front_median": np.median(front) if len(front) else np.nan,
        "core_median": np.median(core) if len(core) else np.nan,
        "front_mean": np.mean(front) if len(front) else np.nan,
        "core_mean": np.mean(core) if len(core) else np.nan,
        "mw_u_pvalue": p
    }

    all_rows.append(row)
    print(row)

# save summary
summary = pd.DataFrame(all_rows)
summary.to_csv(
    f"{OUT_BASE}/REPORTS/front_core_summary_all_samples.csv",
    index=False
)

print("\nDONE ✅")
print("Saved:",
      f"{OUT_BASE}/REPORTS/front_core_summary_all_samples.csv")
