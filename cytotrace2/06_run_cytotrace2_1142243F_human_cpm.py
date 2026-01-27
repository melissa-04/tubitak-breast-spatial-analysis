# ============================================
# STEP13G — Run CytoTRACE2 safely (human, CPM input)
# - detects supported kwargs
# - runs CytoTRACE2
# - extracts scores from returned object OR from output CSV files
# - saves standardized scores CSV to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, glob
import pandas as pd
from inspect import signature
from cytotrace2_py.cytotrace2_py import cytotrace2

ROOT = "/content/drive/MyDrive/TÜBİTAK"
IN_TXT = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RUN_INPUTS/1142243F/cytotrace2_input_1142243F_CPM.txt"

OUT_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RESULTS/1142243F"
CT2_OUTDIR = f"{OUT_DIR}/cytotrace2_results_human"
os.makedirs(CT2_OUTDIR, exist_ok=True)

OUT_CSV = f"{OUT_DIR}/cytotrace2_scores_1142243F_human.csv"

print("Input exists?:", os.path.exists(IN_TXT))

# 1) Build kwargs based on actual signature
params = set(signature(cytotrace2).parameters.keys())
kwargs = {"input_path": IN_TXT, "species": "human"}

if "output_dir" in params:
    kwargs["output_dir"] = CT2_OUTDIR
if "verbose" in params:
    kwargs["verbose"] = True
if "plotting" in params:
    kwargs["plotting"] = False
elif "plotting_enabled" in params:
    kwargs["plotting_enabled"] = False

print("Running with kwargs:", kwargs)

# 2) Run
res = cytotrace2(**kwargs)
print("CytoTRACE2 finished ✅")

# 3) Try to extract scores from return object
score = None
keys_to_try = ["potency_score", "cytotrace2_score", "score", "scores"]

if isinstance(res, dict):
    for k in keys_to_try:
        if k in res:
            score = res[k]
            break
else:
    for k in keys_to_try:
        if hasattr(res, k):
            score = getattr(res, k)
            break

# 4) If not found, search output_dir for a csv containing potency/score
if score is None:
    csvs = glob.glob(os.path.join(CT2_OUTDIR, "**", "*.csv"), recursive=True)
    print("CSV found in output_dir:", len(csvs))
    for p in csvs[:15]:
        print(" -", p)

    if len(csvs) == 0:
        raise ValueError("No CSV found in CytoTRACE2 output_dir; cannot extract scores.")

    # choose best candidate
    priority = [p for p in csvs if ("potency" in p.lower()) or ("score" in p.lower()) or ("cytotrace" in p.lower())]
    pick = priority[0] if len(priority) > 0 else csvs[0]
    print("Using CSV:", pick)

    df = pd.read_csv(pick)
    # if there is a clear score column, keep it; otherwise save whole df
    score_cols = [c for c in df.columns if ("potency" in str(c).lower()) or ("score" in str(c).lower())]
    if len(score_cols) >= 1:
        out = df[[score_cols[0]]].copy()
        out.columns = ["cytotrace2_potency_score"]
    else:
        out = df.copy()
else:
    out = pd.Series(score, name="cytotrace2_potency_score").to_frame()

# 5) Save standardized CSV to Drive
out.index.name = "spot_id"
out.to_csv(OUT_CSV)

print("\nSAVED TO DRIVE ✅ (scores)")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("Preview:")
print(out.head())
# ============================================
# STEP13H — Run CytoTRACE2 (human) + extract scores robustly
# - uses CPM input TXT (genes x spots)
# - forces output_dir to Drive if supported
# - finds potency/score table in outputs
# - saves standardized scores CSV to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, glob
import pandas as pd
from inspect import signature
from cytotrace2_py.cytotrace2_py import cytotrace2

ROOT = "/content/drive/MyDrive/TÜBİTAK"
IN_TXT = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RUN_INPUTS/1142243F/cytotrace2_input_1142243F_CPM.txt"

OUT_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RESULTS/1142243F"
CT2_OUT = f"{OUT_DIR}/ct2_out_human"
os.makedirs(CT2_OUT, exist_ok=True)

OUT_CSV = f"{OUT_DIR}/cytotrace2_scores_1142243F_human.csv"

print("Input exists?:", os.path.exists(IN_TXT))
print("Output dir:", CT2_OUT)

# ---- 1) Build kwargs that THIS version supports
params = set(signature(cytotrace2).parameters.keys())
kwargs = {"input_path": IN_TXT}

# species param name/value can differ; try safest options
if "species" in params:
    kwargs["species"] = "human"

# output dir param name can differ; set if available
for k in ["output_dir", "output_directory", "outdir", "out_dir"]:
    if k in params:
        kwargs[k] = CT2_OUT
        break

# verbose if available
if "verbose" in params:
    kwargs["verbose"] = True

# turn off plotting if supported (name varies)
for k in ["plotting", "plotting_enabled", "plotting_on", "plot"]:
    if k in params:
        kwargs[k] = False
        break

print("Running with kwargs:", kwargs)

# ---- 2) Run CytoTRACE2
res = cytotrace2(**kwargs)
print("\nCytoTRACE2 finished ✅")

# ---- 3) Find output tables (csv/tsv/txt) in output_dir
tables = []
for ext in ["*.csv", "*.tsv", "*.txt"]:
    tables += glob.glob(os.path.join(CT2_OUT, "**", ext), recursive=True)

print("Tables found in output_dir:", len(tables))
for p in tables[:20]:
    print(" -", p)

# If no tables found, fallback to default local output folder
if len(tables) == 0:
    local_dirs = glob.glob("/content/cytotrace2_results*")
    print("\nNo tables in Drive output_dir. Local dirs:", local_dirs)
    for d in local_dirs:
        for ext in ["*.csv", "*.tsv", "*.txt"]:
            tables += glob.glob(os.path.join(d, "**", ext), recursive=True)

    print("Tables found in local dirs:", len(tables))
    for p in tables[:20]:
        print(" -", p)

if len(tables) == 0:
    raise ValueError("No output table found. Cannot extract scores.")

# ---- 4) Pick best candidate table (potency/score keywords)
priority = [p for p in tables if any(k in p.lower() for k in ["potency", "score", "cytotrace"])]
pick = priority[0] if len(priority) > 0 else tables[0]
print("\nUsing table:", pick)

# Read table (auto sep)
if pick.endswith(".csv"):
    df = pd.read_csv(pick)
else:
    df = pd.read_csv(pick, sep="\t")

# ---- 5) Choose a score column
score_cols = [c for c in df.columns if any(k in str(c).lower() for k in ["potency", "cytotrace", "score"])]
if len(score_cols) == 0:
    # last resort: keep full table
    out = df.copy()
else:
    out = df[[score_cols[0]]].copy()
    out.columns = ["cytotrace2_potency_score"]

# If there is an ID column, set it as index
id_cols = [c for c in df.columns if str(c).lower() in ["cell", "cell_id", "barcode", "spot", "spot_id", "id", "unnamed: 0"]]
if len(id_cols) > 0:
    out.index = df[id_cols[0]].astype(str)
out.index.name = "spot_id"

# ---- 6) Save standardized CSV to Drive
out.to_csv(OUT_CSV)
print("\nSAVED TO DRIVE ✅ (scores)")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("Preview:")
print(out.head())
