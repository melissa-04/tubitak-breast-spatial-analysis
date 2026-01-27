# ============================================
# STEP10D — Spatial deconvolution (ONE sample test)
# - loads trained RNA model
# - runs SpatialStereoscope on 1142243F
# - saves proportions (CSV) + updated h5ad to Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import re
import pandas as pd
import scanpy as sc
from scvi.external import RNAStereoscope, SpatialStereoscope

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Paths (inputs)
MODEL_DIR   = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/MODELS/rna_model"
SCRNA_TRAIN = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/INPUTS/scrna_train_downsampled.h5ad"
SPATIAL_IN  = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/INPUTS/spatial_1142243F_counts_intersect.h5ad"

# # Paths (outputs)
OUT_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/RESULTS/1142243F"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_CSV  = f"{OUT_DIR}/stereoscope_proportions_1142243F.csv"
OUT_H5AD = f"{OUT_DIR}/spatial_1142243F_stereoscope_result.h5ad"

print("MODEL_DIR exists?:", os.path.exists(MODEL_DIR))
print("SCRNA_TRAIN exists?:", os.path.exists(SCRNA_TRAIN))
print("SPATIAL_IN exists?:", os.path.exists(SPATIAL_IN))

# # 1) Load scRNA training data (needed to load the model registry)
adata_sc = sc.read_h5ad(SCRNA_TRAIN)

# # counts layer ekle (Stereoscope setup için)
adata_sc.layers["counts"] = adata_sc.X.copy()

# # setup (labels_key = cell_type)
RNAStereoscope.setup_anndata(adata_sc, layer="counts", labels_key="cell_type")

# # 2) Load trained RNA model from Drive
rna_model = RNAStereoscope.load(MODEL_DIR, adata=adata_sc)
print("\nRNA model loaded ✅")

# # 3) Load spatial sample
adata_sp = sc.read_h5ad(SPATIAL_IN)

# # counts layer ekle (spatial model setup için)
adata_sp.layers["counts"] = adata_sp.X.copy()

# # setup for spatial (no labels)
SpatialStereoscope.setup_anndata(adata_sp, layer="counts")

# # 4) Build spatial model from RNA model and train
sp_model = SpatialStereoscope.from_rna_model(adata_sp, rna_model)
sp_model.train(max_epochs=200)

print("\nSpatial model trained ✅")

# # 5) Get proportions (spots x cell_types)
props = sp_model.get_proportions()
# props genelde DataFrame gelir; değilse DataFrame'e çevir
if not isinstance(props, pd.DataFrame):
    props = pd.DataFrame(props)

print("\nProportions shape:", props.shape)
print("Columns (cell types):", list(props.columns))

# # 6) Save proportions to Drive (CSV)
props.to_csv(OUT_CSV, index=True)
print("\nSAVED TO DRIVE ✅ (proportions CSV)")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))

# # 7) Add proportions into adata_sp (obs + obsm) and save h5ad
# # obs kolon isimlerini güvenli hale getir
def safe_name(x):
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    return x.strip("_")

adata_sp.obsm["stereoscope_proportions"] = props.values
adata_sp.uns["stereoscope_cell_types"] = list(props.columns)

for col in props.columns:
    adata_sp.obs[f"stereo_{safe_name(col)}"] = props[col].values

adata_sp.write_h5ad(OUT_H5AD)

print("\nSAVED TO DRIVE ✅ (spatial h5ad with results)")
print("Path:", OUT_H5AD)
print("Saved OK?:", os.path.exists(OUT_H5AD))

print("\nDONE ✅")
