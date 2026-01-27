# ============================================
# STEP10C (FIX) — Train RNA Stereoscope model
# Uses scvi.external.RNAStereoscope (correct class)
# Saves model + training history to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import scanpy as sc
from scvi.external import RNAStereoscope  # # doğru import

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # INPUT: bir önceki adımda kaydettiğin downsample scRNA
SCRNA_TRAIN = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/INPUTS/scrna_train_downsampled.h5ad"

# # OUTPUT: model + rapor
MODEL_DIR  = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/MODELS/rna_model"
REPORT_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/REPORTS"
os.makedirs(MODEL_DIR, exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)

# # 1) Load training data
adata = sc.read_h5ad(SCRNA_TRAIN)
print("Training scRNA shape:", adata.shape)

# # 2) Stereoscope counts layer ister -> X'i counts layer'a koy
adata.layers["counts"] = adata.X.copy()

# # 3) Setup + Train (labels_key = cell_type)
RNAStereoscope.setup_anndata(adata, layer="counts", labels_key="cell_type")
sc_model = RNAStereoscope(adata)
sc_model.train(max_epochs=100)

# # 4) Save model to Drive
sc_model.save(MODEL_DIR, overwrite=True)
print("\nSAVED TO DRIVE ✅ (RNA model)")
print("Path:", MODEL_DIR)
print("Model dir exists?:", os.path.exists(MODEL_DIR))

# # 5) Save training history to Drive
hist_df = pd.DataFrame(sc_model.history)
hist_path = f"{REPORT_DIR}/rna_model_training_history.csv"
hist_df.to_csv(hist_path, index=False)

print("\nSAVED TO DRIVE ✅ (training history)")
print("Path:", hist_path)
print("Saved OK?:", os.path.exists(hist_path))

print("\nDONE: RNA model training completed.")
