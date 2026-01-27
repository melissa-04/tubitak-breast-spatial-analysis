# ============================================
# STEP08 — Normalize QC2_filtered: CPM -> log2
# Saves outputs to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import scanpy as sc

# # Proje kökü (Drive)
PROJECT_ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Girdi: QC2_filtered h5ad'ler
QC2_DIR = f"{PROJECT_ROOT}/outputs/STEP07B_H5AD_QC_FILTERED"

# # Çıktı: normalizasyon sonucu (Drive'a yazılacak)
OUT_DIR = f"{PROJECT_ROOT}/outputs/STEP08_NORMALIZATION/QC2"
os.makedirs(OUT_DIR, exist_ok=True)

print("QC2 input dir :", QC2_DIR)
print("QC2 output dir:", OUT_DIR)

# # QC2 dosyalarını bul (alt klasörler dahil)
qc2_paths = []
for root, dirs, files in os.walk(QC2_DIR):
    for f in files:
        if f.endswith(".h5ad") and "QC2_filtered" in f:
            qc2_paths.append(os.path.join(root, f))

qc2_paths = sorted(qc2_paths)
print("\nFound QC2 h5ad files:", len(qc2_paths))
for p in qc2_paths:
    print(" -", p)

def normalize_to_log2cpm(in_path, out_path, target_sum=1e6):
    # # Oku
    adata = sc.read_h5ad(in_path)

    # # Ham sayımı sakla
    adata.layers["counts"] = adata.X.copy()

    # # CPM normalize (spot bazlı)
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

    # # log2(1+x)
    sc.pp.log1p(adata, base=2)

    # # Normalize matrisi ayrı katmanda sakla
    adata.layers["log2_cpm"] = adata.X.copy()

    # # Ne yaptığımızı kaydet
    adata.uns["normalization_step08"] = {
        "method": "CPM + log2(1+x)",
        "target_sum": float(target_sum),
        "raw_layer": "counts",
        "norm_layer": "log2_cpm"
    }

    # # Drive'a kaydet
    adata.write_h5ad(out_path)

    # # Kaydedildi mi kontrol
    print("\nSaved:", out_path)
    print("Saved OK?:", os.path.exists(out_path))
    print("Shape:", adata.shape)
    print("Layers:", list(adata.layers.keys()))

# # Tüm QC2 dosyalarını normalize et
for in_path in qc2_paths:
    base = os.path.basename(in_path).replace(".h5ad", "")
    out_path = f"{OUT_DIR}/{base}__log2cpm.h5ad"
    normalize_to_log2cpm(in_path, out_path)

print("\nDONE: QC2 normalization finished.")
