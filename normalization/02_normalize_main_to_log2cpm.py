# ============================================
# STEP08B — Normalize MAIN: CPM -> log2
# Saves outputs to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import scanpy as sc
import numpy as np

# # Proje kökü
PROJECT_ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Girdi: MAIN h5ad'ler (artefact yok, in_tissue)
MAIN_DIR = f"{PROJECT_ROOT}/outputs/STEP07_H5AD_MAIN"

# # Çıktı: MAIN normalizasyon sonucu
OUT_DIR = f"{PROJECT_ROOT}/outputs/STEP08_NORMALIZATION/MAIN"
os.makedirs(OUT_DIR, exist_ok=True)

print("MAIN input dir :", MAIN_DIR)
print("MAIN output dir:", OUT_DIR)

# # MAIN dosyalarını bul (alt klasörler dahil)
main_paths = []
for root, dirs, files in os.walk(MAIN_DIR):
    for f in files:
        if f.endswith(".h5ad") and "MAIN_keep_in_tissue_no_artefact" in f:
            main_paths.append(os.path.join(root, f))

main_paths = sorted(main_paths)
print("\nFound MAIN h5ad files:", len(main_paths))
for p in main_paths:
    print(" -", p)

def normalize_to_log2cpm(in_path, out_path, target_sum=1e6):
    # # Oku
    adata = sc.read_h5ad(in_path)

    # # Ham sayımı sakla
    adata.layers["counts"] = adata.X.copy()

    # # Küçük kontrol: counts dtype (ileride count-based yöntemler için önemli olabilir)
    try:
        dt = adata.layers["counts"].dtype
    except Exception:
        dt = "unknown"
    print("\nInput:", os.path.basename(in_path))
    print("counts dtype:", dt)

    # # CPM normalize
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

    # # log2(1+x)
    sc.pp.log1p(adata, base=2)

    # # Normalize matrisi sakla
    adata.layers["log2_cpm"] = adata.X.copy()

    # # Metaya yaz
    adata.uns["normalization_step08"] = {
        "method": "CPM + log2(1+x)",
        "target_sum": float(target_sum),
        "raw_layer": "counts",
        "norm_layer": "log2_cpm",
        "input_type": "MAIN_keep_in_tissue_no_artefact"
    }

    # # Drive'a kaydet
    adata.write_h5ad(out_path)

    # # Kaydedildi mi?
    print("Saved:", out_path)
    print("Saved OK?:", os.path.exists(out_path))
    print("Shape:", adata.shape)
    print("Layers:", list(adata.layers.keys()))

# # Hepsini çalıştır
for in_path in main_paths:
    base = os.path.basename(in_path).replace(".h5ad", "")
    out_path = f"{OUT_DIR}/{base}__log2cpm.h5ad"
    normalize_to_log2cpm(in_path, out_path)

print("\nDONE: MAIN normalization finished.")
