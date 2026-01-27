# ============================================
# STEP13C — Build CytoTRACE2 input TXT (genes x spots) for 1142243F
# Uses log2_cpm layer from cytotrace2-ready h5ad
# Saves TXT to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Input: cytotrace2-ready strict set (has layers['counts'] and layers['log2_cpm'])
IN_H5AD = f"{ROOT}/outputs/STEP13_CYTOTRACE2/INPUTS/1142243F/spatial_1142243F__strict__cytotrace2_ready.h5ad"

# # Output TXT
OUT_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/RUN_INPUTS/1142243F"
os.makedirs(OUT_DIR, exist_ok=True)
OUT_TXT = f"{OUT_DIR}/cytotrace2_input_1142243F.txt"

adata = sc.read_h5ad(IN_H5AD)

# # Use log2_cpm as expression matrix for CytoTRACE2
X = adata.layers["log2_cpm"]
if issparse(X):
    X = X.toarray()

# # We need genes as rows, spots as columns -> transpose
# adata: (spots x genes)  ->  (genes x spots)
expr = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)

# # Save tab-delimited
expr.to_csv(OUT_TXT, sep="\t")

print("SAVED TO DRIVE ✅ (CytoTRACE2 input)")
print("Path:", OUT_TXT)
print("Saved OK?:", os.path.exists(OUT_TXT))
print("Matrix shape (genes x spots):", expr.shape)
