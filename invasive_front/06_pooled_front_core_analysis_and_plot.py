# ============================================
# STEP16C — Pooled front vs core analysis + violin plot
# - loads per-sample spot tables: candidate_front_core_scored_<SAMPLE>.csv
# - pools all front spots and all core spots
# - runs MWU test
# - saves stats + plot to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = "/content/drive/MyDrive/TÜBİTAK"
RES_BASE = f"{ROOT}/outputs/STEP16_CANDIDATE_CYTOTRACE2/RESULTS"
REP_DIR = f"{ROOT}/outputs/STEP16_CANDIDATE_CYTOTRACE2/REPORTS"
FIG_DIR = f"{ROOT}/outputs/STEP16_CANDIDATE_CYTOTRACE2/FIGURES"
os.makedirs(REP_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
KEY = "cytotrace2_potency_score"

dfs = []
for s in samples:
    fp = f"{RES_BASE}/{s}/candidate_front_core_scored_{s}.csv"
    df = pd.read_csv(fp, index_col=0)
    df["sample"] = s
    dfs.append(df)

data = pd.concat(dfs, axis=0)

front = data[data["region"]=="front"][KEY].astype(float).values
core  = data[data["region"]=="core"][KEY].astype(float).values

u, p = mannwhitneyu(front, core, alternative="two-sided")

stats = pd.DataFrame([{
    "n_front": len(front),
    "n_core": len(core),
    "front_median": float(np.median(front)),
    "core_median": float(np.median(core)),
    "delta_median": float(np.median(front) - np.median(core)),
    "front_mean": float(np.mean(front)),
    "core_mean": float(np.mean(core)),
    "delta_mean": float(np.mean(front) - np.mean(core)),
    "mwu_pvalue": float(p)
}])

out_stats = f"{REP_DIR}/pooled_front_core_stats.csv"
stats.to_csv(out_stats, index=False)

print("SAVED TO DRIVE ✅ (pooled stats)")
print("Path:", out_stats)
print("Saved OK?:", os.path.exists(out_stats))
print(stats)

# Plot (violin style using matplotlib: simple boxplot to keep it clean)
out_png = f"{FIG_DIR}/pooled_front_vs_core_violin.png"

plt.figure()
plt.boxplot([front, core], labels=["Front", "Core"])
plt.ylabel("CytoTRACE2 potency score")
plt.title(f"Pooled Front vs Core (MWU p={p:.2e})")
plt.savefig(out_png, dpi=220, bbox_inches="tight")
plt.close()

print("\nSAVED TO DRIVE ✅ (figure)")
print("Path:", out_png)
print("Saved OK?:", os.path.exists(out_png))
