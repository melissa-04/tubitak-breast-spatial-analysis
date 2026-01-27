# ============================================
# STEP13I — Run CytoTRACE2 on ALL strict samples (CPM, human)
# - builds CPM input (genes x spots)
# - runs cytotrace2 with signature-safe kwargs
# - extracts potency scores from output tables
# - saves per-sample scores CSV + run summary to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from inspect import signature
from cytotrace2_py.cytotrace2_py import cytotrace2

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# Inputs: cytotrace2-ready strict sets (have layers['counts'])
IN_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/INPUTS"

# Outputs
RUN_INPUT_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RUN_INPUTS"
RES_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RESULTS"
REPORT_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/REPORTS"
os.makedirs(RUN_INPUT_BASE, exist_ok=True)
os.makedirs(RES_BASE, exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

# Detect supported kwargs once
params = set(signature(cytotrace2).parameters.keys())

def run_one(sample):
    # Paths
    in_h5ad = f"{IN_BASE}/{sample}/spatial_{sample}__strict__cytotrace2_ready.h5ad"

    run_in_dir = f"{RUN_INPUT_BASE}/{sample}"
    res_dir = f"{RES_BASE}/{sample}"
    ct2_out = f"{res_dir}/ct2_out_human"
    os.makedirs(run_in_dir, exist_ok=True)
    os.makedirs(res_dir, exist_ok=True)
    os.makedirs(ct2_out, exist_ok=True)

    in_txt = f"{run_in_dir}/cytotrace2_input_{sample}_CPM.txt"
    out_csv = f"{res_dir}/cytotrace2_scores_{sample}_human.csv"

    # Skip if already exists
    if os.path.exists(out_csv):
        return {
            "sample": sample,
            "status": "SKIPPED_EXISTS",
            "scores_csv": out_csv,
            "n_spots": None,
            "n_genes": None
        }

    # ---- Build CPM input from counts
    adata = sc.read_h5ad(in_h5ad)
    X = adata.layers["counts"]
    if issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)

    row_sums = X.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    X_cpm = (X / row_sums) * 1e6

    expr = pd.DataFrame(X_cpm.T, index=adata.var_names, columns=adata.obs_names)  # genes x spots
    expr.to_csv(in_txt, sep="\t")

    # ---- Run cytotrace2 with safe kwargs
    kwargs = {"input_path": in_txt}
    if "species" in params:
        kwargs["species"] = "human"
    for k in ["output_dir", "output_directory", "outdir", "out_dir"]:
        if k in params:
            kwargs[k] = ct2_out
            break
    if "verbose" in params:
        kwargs["verbose"] = True
    for k in ["plotting", "plotting_enabled", "plotting_on", "plot"]:
        if k in params:
            kwargs[k] = False
            break

    res = cytotrace2(**kwargs)

    # ---- Find output tables and extract score
    tables = []
    for ext in ["*.csv", "*.tsv", "*.txt"]:
        tables += glob.glob(os.path.join(ct2_out, "**", ext), recursive=True)

    if len(tables) == 0:
        # fallback: local default dir
        for d in glob.glob("/content/cytotrace2_results*"):
            for ext in ["*.csv", "*.tsv", "*.txt"]:
                tables += glob.glob(os.path.join(d, "**", ext), recursive=True)

    if len(tables) == 0:
        return {
            "sample": sample,
            "status": "FAILED_NO_TABLE",
            "scores_csv": "",
            "n_spots": adata.n_obs,
            "n_genes": adata.n_vars
        }

    priority = [p for p in tables if any(k in p.lower() for k in ["potency", "score", "cytotrace"])]
    pick = priority[0] if len(priority) > 0 else tables[0]

    if pick.endswith(".csv"):
        df = pd.read_csv(pick)
    else:
        df = pd.read_csv(pick, sep="\t")

    score_cols = [c for c in df.columns if any(k in str(c).lower() for k in ["potency", "cytotrace", "score"])]
    if len(score_cols) >= 1:
        out = df[[score_cols[0]]].copy()
        out.columns = ["cytotrace2_potency_score"]
    else:
        # last resort
        out = df.copy()

    id_cols = [c for c in df.columns if str(c).lower() in ["cell", "cell_id", "barcode", "spot", "spot_id", "id", "unnamed: 0"]]
    if len(id_cols) > 0 and "cytotrace2_potency_score" in out.columns:
        out.index = df[id_cols[0]].astype(str)

    out.index.name = "spot_id"
    out.to_csv(out_csv)

    return {
        "sample": sample,
        "status": "OK",
        "scores_csv": out_csv,
        "n_spots": adata.n_obs,
        "n_genes": adata.n_vars
    }

# Run all
summary = []
for s in samples:
    print(f"\n=== {s} ===")
    info = run_one(s)
    summary.append(info)
    print("Status:", info["status"])
    print("Scores CSV:", info["scores_csv"])
    if info["status"] == "OK":
        print("Saved OK?:", os.path.exists(info["scores_csv"]))

# Save summary report
summary_df = pd.DataFrame(summary)
sum_path = f"{REPORT_DIR}/cytotrace2_all_samples_summary.csv"
summary_df.to_csv(sum_path, index=False)

print("\nSAVED TO DRIVE ✅ (summary)")
print("Path:", sum_path)
print("Saved OK?:", os.path.exists(sum_path))
print(summary_df)
