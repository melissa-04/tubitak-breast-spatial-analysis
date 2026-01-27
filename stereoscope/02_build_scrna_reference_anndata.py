# ============================================
# STEP09C (FIX) — Build scRNA reference AnnData (counts)
# Fix: transpose mtx if it is genes x cells
# Output SAVED to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix

# # Proje kökü
ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Extract edilen scRNA klasörü
SCRNA_DIR = f"{ROOT}/reference/GSE176078/scrna/Wu_etal_2021_BRCA_scRNASeq"

# # Girdi dosyaları
MTX  = f"{SCRNA_DIR}/count_matrix_sparse.mtx"
GENE = f"{SCRNA_DIR}/count_matrix_genes.tsv"
BC   = f"{SCRNA_DIR}/count_matrix_barcodes.tsv"
META = f"{SCRNA_DIR}/metadata.csv"

# # 1) Dosyaları oku
X = mmread(MTX).tocsr()  # # sparse counts
genes = pd.read_csv(GENE, sep="\t", header=None)[0].astype(str).values
barcodes = pd.read_csv(BC, sep="\t", header=None)[0].astype(str).values
meta = pd.read_csv(META)

print("Loaded:")
print(" - raw X shape:", X.shape)
print(" - genes:", len(genes))
print(" - barcodes:", len(barcodes))
print(" - meta rows/cols:", meta.shape)

# # 2) Eğer X (genes x cells) ise transpose yap -> (cells x genes)
# # Kural: X.shape == (len(barcodes), len(genes)) olmalı
if X.shape == (len(genes), len(barcodes)):
    print("\nTransposing X: (genes x cells) -> (cells x genes)")
    X = X.T.tocsr()
elif X.shape == (len(barcodes), len(genes)):
    print("\nX already in (cells x genes) format ✅")
else:
    raise ValueError(f"Unexpected X shape: {X.shape}. genes={len(genes)}, barcodes={len(barcodes)}")

print(" - fixed X shape:", X.shape)

# # 3) AnnData oluştur (obs=cells, var=genes)
adata = sc.AnnData(X=csr_matrix(X))
adata.obs_names = barcodes
adata.var_names = genes

# # 4) metadata'da barcode/cell id kolonu var mı? (adayları yazdır)
candidates = [c for c in meta.columns if str(c).lower() in
              ["barcode", "barcodes", "cell", "cell_id", "cellid", "id", "unnamed: 0"]]
print("\nMetadata columns:", list(meta.columns))
print("Barcode-like candidate columns:", candidates)

# # 5) metadata index'i barcode yap ve hizala
if len(candidates) >= 1:
    key = candidates[0]
    meta[key] = meta[key].astype(str)
    meta = meta.set_index(key)
else:
    meta.index = meta.index.astype(str)

meta_aligned = meta.reindex(adata.obs_names)
adata.obs = meta_aligned

# # 6) Counts'u int yap (Stereoscope için güvenli)
adata.X = adata.X.tocsr()
adata.X.data = np.rint(adata.X.data).astype(np.int64)

# # 7) Kontroller
print("\nSanity checks:")
print(" - adata shape:", adata.shape)
print(" - X dtype:", adata.X.data.dtype)
print(" - total counts:", int(adata.X.sum()))
print(" - cells with ALL metadata missing:", int(adata.obs.isna().all(axis=1).sum()))

# # 8) Label kolon adaylarını yazdır (sonraki adımda seçeceğiz)
cols = list(adata.obs.columns)
type_like = [c for c in cols if any(k in str(c).lower() for k in ["type", "annot", "cluster", "subtype", "label"])]
print("\nPossible label columns (candidates):")
for c in type_like[:40]:
    print(" -", c)

# # 9) Drive'a KAYDET
OUT_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE"
os.makedirs(OUT_DIR, exist_ok=True)

out_path = f"{OUT_DIR}/GSE176078_scRNA_reference_counts.h5ad"
adata.write_h5ad(out_path)

print("\nSAVED TO DRIVE ✅")
print("Path:", out_path)
print("Saved OK?:", os.path.exists(out_path))
