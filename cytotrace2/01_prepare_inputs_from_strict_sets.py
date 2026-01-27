# ============================================
# STEP13A — Prepare CytoTRACE2-ready inputs from STRICT sets
# - ensures counts + log2_cpm layers exist
# - saves per-sample cytotrace2-ready h5ad to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Strict sets
STRICT_BASE = f"{ROOT}/outputs/STEP12_STRICT_CANCER_RICH"

# # Output
OUT_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/INPUTS"
REPORT_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/REPORTS"
os.makedirs(OUT_BASE, exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

rows = []

for s in samples:
    in_h5ad = f"{STRICT_BASE}/{s}/spatial_{s}__strict_cancer_rich.h5ad"
    adata = sc.read_h5ad(in_h5ad)

    # # 1) counts layer: ensure integer counts
    adata.layers["counts"] = adata.X.copy()
    if hasattr(adata.layers["counts"], "data"):
        adata.layers["counts"].data = np.rint(adata.layers["counts"].data).astype(np.int64)
    else:
        adata.layers["counts"] = np.rint(adata.layers["counts"]).astype(np.int64)

    # # 2) build log2 CPM in X, then store as log2_cpm
    sc.pp.normalize_total(adata, target_sum=1e6, inplace=True)  # CPM
    sc.pp.log1p(adata, base=2)                                  # log2(1+x)
    adata.layers["log2_cpm"] = adata.X.copy()

    # # 3) Save cytotrace2-ready file
    out_dir = f"{OUT_BASE}/{s}"
    os.makedirs(out_dir, exist_ok=True)
    out_path = f"{out_dir}/spatial_{s}__strict__cytotrace2_ready.h5ad"
    adata.write_h5ad(out_path)

    rows.append({
        "sample": s,
        "n_spots": adata.n_obs,
        "n_genes": adata.n_vars,
        "out_path": out_path
    })

    print(f"\n{s}")
    print(" - spots:", adata.n_obs, "genes:", adata.n_vars)
    print(" - layers:", list(adata.layers.keys()))
    print(" - saved OK?:", os.path.exists(out_path))

# # Save summary
df = pd.DataFrame(rows)
summary_path = f"{REPORT_DIR}/cytotrace2_input_summary.csv"
df.to_csv(summary_path, index=False)

print("\nSAVED TO DRIVE ✅ (summary)")
print("Path:", summary_path)
print("Saved OK?:", os.path.exists(summary_path))
print(df)
