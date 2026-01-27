# ============================================
# STEP13E — Fix environment for CytoTRACE2 (clean install) + verify
# Saves package versions to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, importlib.metadata as md

# # 1) Install a consistent stack (NumPy>=2 required by several deps)
!pip -q install -U \
  "numpy>=2.0" \
  "scipy>=1.12" \
  "pandas>=2.0" \
  "h5py>=3.10" \
  "anndata>=0.10.8" \
  "scanpy>=1.10" \
  "scikit-learn>=1.4" \
  "umap-learn>=0.5.6" \
  "cytotrace2-py"

# # 2) Import tests
import numpy as np
import scipy
import pandas as pd
import scanpy as sc
from cytotrace2_py.cytotrace2_py import cytotrace2

print("OK ✅ Imports successful")
print("numpy:", np.__version__)
print("scipy:", scipy.__version__)
print("pandas:", pd.__version__)
print("scanpy:", md.version("scanpy"))
print("cytotrace2-py:", md.version("cytotrace2-py"))

# # 3) Save versions to Drive (reproducibility)
ROOT = "/content/drive/MyDrive/TÜBİTAK"
REP_DIR = f"{ROOT}/outputs/STEP13_CYTOTRACE2/REPORTS"
os.makedirs(REP_DIR, exist_ok=True)
out_path = f"{REP_DIR}/env_versions_cytotrace2.txt"

with open(out_path, "w") as w:
    w.write(f"numpy={np.__version__}\n")
    w.write(f"scipy={scipy.__version__}\n")
    w.write(f"pandas={pd.__version__}\n")
    w.write(f"scanpy={md.version('scanpy')}\n")
    w.write(f"cytotrace2-py={md.version('cytotrace2-py')}\n")

print("\nSAVED TO DRIVE ✅")
print("Path:", out_path)
print("Saved OK?:", os.path.exists(out_path))
