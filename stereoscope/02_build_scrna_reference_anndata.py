# ============================================
# STEP09C — Build scRNA reference AnnData (counts)
# Output is SAVED to Google Drive at the end
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

# # Extract edilen scRNA klasörü (senin çıktına göre)
SCRNA_DIR = f"{ROOT}/reference/GSE176078/scrna/Wu_etal_2021_BRCA_scRNASeq"

# # Girdi dosyaları
MTX  = f"{SCRNA_DIR}/count_matrix_sparse.mtx"
GENE = f"{SCRNA_DIR}/count_matrix_genes.tsv"
BC   = f"{SCRNA_DIR}/count_matrix_barcodes.tsv"
META = f"{SCRNA_DIR}/metadata.csv"

print("SCRNA_DIR:", SCRNA_DIR)
print("Exists MTX? :", os.path.exists(MTX))
print("Exists GENE?:", os.path.exists(GENE))
print("Exists BC?  :", os.path.exists(BC))
print("Exists META?:", os.path.exists(META))

# # 1) Count matrix + gene + barcode + metadata oku
X = mmread(MTX).tocsr()  # # sparse counts
genes = pd.read_csv(GENE, sep="\t", header=None)[0].astype(str).values
barcodes = pd.read_csv(BC, sep="\t", header=None)[0].astype(str).values
meta = pd.read_csv(META)

print("\nLoaded:")
print(" - X shape:", X.shape)
print(" - genes   :", len(genes))
print(" - barcodes:", len(barcodes))
print(" - meta rows/cols:", meta.shape)

# # 2) AnnData oluştur
adata = sc.AnnData(X=csr_matrix(X))
adata.var_names = genes
adata.obs_names = barcodes

# # 3) metadata'da barcode/cell id hangi kolonda? (adayları yazdır)
candidates = [c for c in meta.columns if str(c).lower() in
              ["barcode", "barcodes", "cell", "cell_id", "cellid", "id", "unnamed: 0"]]
print("\nMetadata columns (first 20):", list(meta.columns)[:20])
print("Barcode-like candidate columns:", candidates)

# # 4) metadata index'i barcode yap ve adata.obs'a hizala
if len(candidates) >= 1:
    key = candidates[0]  # # ilk adayı kullanıyoruz (gerekirse sonra değiştiririz)
    meta[key] = meta[key].astype(str)
    meta = meta.set_index(key)
else:
    # # aday yoksa: metadata index zaten barcode olabilir
    meta.index = meta.index.astype(str)

meta_aligned = meta.reindex(adata.obs_names)
adata.obs = meta_aligned

# # 5) Counts'u int yap (Stereoscope için güvenli)
adata.X = adata.X.tocsr()
adata.X.data = np.rint(adata.X.data).astype(np.int64)

# # 6) Kontroller
print("\nSanity checks:")
print(" - adata shape:", adata.shape)
print(" - X dtype:", adata.X.data.dtype)
print(" - total counts:", int(adata.X.sum()))
print(" - cells with ALL metadata missing:", int(adata.obs.isna().all(axis=1).sum()))

# # 7) Cell-type benzeri kolonları öner (sonraki adımda doğru kolon seçip standartlaştıracağız)
cols = list(adata.obs.columns)
type_like = [c for c in cols if any(k in str(c).lower() for k in ["type", "annot", "cluster", "subtype", "label"])]
print("\nPossible label columns (candidates):")
for c in type_like[:40]:
    print(" -", c)

# # 8) Drive'a KAYDET (nereye kaydettiğini ekrana yaz)
OUT_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE"
os.makedirs(OUT_DIR, exist_ok=True)

out_path = f"{OUT_DIR}/GSE176078_scRNA_reference_counts.h5ad"
adata.write_h5ad(out_path)

print("\nSAVED TO DRIVE ✅")
print("Path:", out_path)
print("Saved OK?:", os.path.exists(out_path))
