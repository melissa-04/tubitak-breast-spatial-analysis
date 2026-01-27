# ============================================
# STEP16B — Add FDR (BH) + effect sizes to per-sample stats
# Saves updated summary CSV to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd

ROOT = "/content/drive/MyDrive/TÜBİTAK"
IN_CSV = f"{ROOT}/outputs/STEP16_CANDIDATE_CYTOTRACE2/REPORTS/front_core_stats_candidate_per_sample.csv"
OUT_CSV = f"{ROOT}/outputs/STEP16_CANDIDATE_CYTOTRACE2/REPORTS/front_core_stats_candidate_with_fdr.csv"

df = pd.read_csv(IN_CSV)

# Effect sizes
df["delta_median"] = df["front_median"] - df["core_median"]
df["delta_mean"]   = df["front_mean"]   - df["core_mean"]

# BH-FDR
p = df["mw_u_pvalue"].astype(float).values
m = len(p)
order = np.argsort(p)
rank = np.empty(m, dtype=int)
rank[order] = np.arange(1, m+1)

fdr = (p * m) / rank
# enforce monotonicity (BH)
fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
fdr_final = np.empty(m)
fdr_final[order] = np.minimum(fdr_sorted, 1.0)

df["fdr_bh"] = fdr_final

# Save
df.to_csv(OUT_CSV, index=False)

print("SAVED TO DRIVE ✅")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print(df[["sample","mw_u_pvalue","fdr_bh","delta_median","delta_mean"]])
