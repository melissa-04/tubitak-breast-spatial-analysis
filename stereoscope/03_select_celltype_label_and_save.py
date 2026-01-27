# ============================================
# STEP09D — Choose cell type label column + save
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

IN_PATH = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE/GSE176078_scRNA_reference_counts.h5ad"
OUT_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE"
os.makedirs(OUT_DIR, exist_ok=True)

# # 1) Load reference
adata = sc.read_h5ad(IN_PATH)

# # 2) Choose label column (recommended: celltype_major)
LABEL_COL = "celltype_major"

# # 3) Create standard label column for downstream tools
adata.obs["cell_type"] = adata.obs[LABEL_COL].astype(str)

# # 4) Quick summary
print("Label column used:", LABEL_COL)
print("Unique cell_type count:", adata.obs["cell_type"].nunique())
print("Top 10 labels:")
print(adata.obs["cell_type"].value_counts().head(10))

# # 5) Save to Drive
OUT_PATH = f"{OUT_DIR}/GSE176078_scRNA_reference_counts__label_{LABEL_COL}.h5ad"
adata.write_h5ad(OUT_PATH)

print("\nSAVED TO DRIVE ✅")
print("Path:", OUT_PATH)
print("Saved OK?:", os.path.exists(OUT_PATH))
