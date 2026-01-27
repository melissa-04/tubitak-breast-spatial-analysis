# ============================================
# STEP13D — Run CytoTRACE2 on 1142243F
# - reads TXT (genes x spots)
# - runs cytotrace2
# - saves spot-level scores to Drive (CSV)
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import pickle

from cytotrace2_py.cytotrace2_py import cytotrace2

ROOT = "/content/drive/MyDrive/TÜBİTAK"

IN_TXT = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RUN_INPUTS/1142243F/cytotrace2_input_1142243F.txt"

OUT_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RESULTS/1142243F"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_CSV = f"{OUT_DIR}/cytotrace2_scores_1142243F.csv"
OUT_PKL = f"{OUT_DIR}/cytotrace2_output_full_1142243F.pkl"

print("Input exists?:", os.path.exists(IN_TXT))
print("Running CytoTRACE2...")

# # Run (returns results object)
res = cytotrace2(input_path=IN_TXT)

print("CytoTRACE2 finished ✅")

# # Save full object (for debugging / reproducibility)
with open(OUT_PKL, "wb") as f:
    pickle.dump(res, f)
print("Saved full output ✅:", OUT_PKL, "OK?", os.path.exists(OUT_PKL))

# # Extract spot-level potency score
# Different versions may store results differently, so we handle common patterns.
score_series = None

# Case 1: dict-like
if isinstance(res, dict):
    # try common keys
    for k in ["potency_score", "cytotrace2_score", "score", "scores"]:
        if k in res:
            score_series = res[k]
            break

# Case 2: attribute-like
if score_series is None:
    for attr in ["potency_score", "cytotrace2_score", "score", "scores"]:
        if hasattr(res, attr):
            score_series = getattr(res, attr)
            break

# If still None, show keys for troubleshooting
if score_series is None:
    print("Could not find score directly. Type:", type(res))
    if isinstance(res, dict):
        print("Keys:", list(res.keys()))
    raise ValueError("CytoTRACE2 output format unexpected; cannot extract scores.")

# Convert to DataFrame
scores = pd.Series(score_series)
scores_df = scores.to_frame(name="cytotrace2_potency_score")
scores_df.index.name = "spot_id"
scores_df.to_csv(OUT_CSV)

print("\nSAVED TO DRIVE ✅ (scores)")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("Preview:")
print(scores_df.head())
