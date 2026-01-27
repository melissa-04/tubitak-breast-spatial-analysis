# ============================================
# STEP08 — Normalize QC2_filtered (CPM -> log2)
# Output files are SAVED to Google Drive
# ============================================

import os
import scanpy as sc

# -------------------------------
# 1) Drive paths (INPUT / OUTPUT)
# -------------------------------
# # PROJECT_ROOT: Google Drive içindeki ana proje klasörün
PROJECT_ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # QC2_DIR: Normalizasyona girecek temiz (QC2_filtered) h5ad dosyalarının bulunduğu yer
QC2_DIR = f"{PROJECT_ROOT}/outputs/STEP07B_H5AD_QC_FILTERED"

# # OUT_DIR: Normalizasyon çıktılarının kaydedileceği klasör (Drive üzerinde)
OUT_DIR = f"{PROJECT_ROOT}/outputs/STEP08_NORMALIZATION/QC2"
os.makedirs(OUT_DIR, exist_ok=True)  # # yoksa oluştur

print("QC2 input dir :", QC2_DIR)
print("QC2 output dir:", OUT_DIR)

# ----------------------------------------
# 2) QC2 h5ad dosyalarını bul (recursive)
# ----------------------------------------
qc2_paths = []
for root, dirs, files in os.walk(QC2_DIR):
    for f in files:
        # # sadece .h5ad ve dosya isminde QC2_filtered geçenleri al
        if f.endswith(".h5ad") and "QC2_filtered" in f:
            qc2_paths.append(os.path.join(root, f))

qc2_paths = sorted(qc2_paths)

print("\nFound QC2 h5ad files:", len(qc2_paths))
for p in qc2_paths:
    print(" -", p)

# ----------------------------------------------------
# 3) Normalizasyon fonksiyonu (counts'u koru!)
# ----------------------------------------------------
def normalize_to_log2cpm(in_path, out_path, target_sum=1e6):
    """
    # Amaç:
    # - ham sayımları kaybetme (layers['counts'])
    # - CPM normalize et (target_sum)
    # - log2(1 + CPM) dönüşümü yap
    # - sonucu layers['log2_cpm'] içine yaz
    # - .h5ad olarak DRIVE'a kaydet
    """

    # # 3.1 Dosyayı oku
    adata = sc.read_h5ad(in_path)

    # # 3.2 Ham sayımı sakla (dokunmadan)
    adata.layers["counts"] = adata.X.copy()

    # # 3.3 CPM normalize (spot başına toplamı target_sum yap)
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)

    # # 3.4 log2(1 + x) dönüşümü
    sc.pp.log1p(adata, base=2)

    # # 3.5 Normalized matrisi ayrı katmanda sakla
    adata.layers["log2_cpm"] = adata.X.copy()

    # # 3.6 Yapılan işlemi kayda geçir (reproducibility)
    adata.uns["normalization_step08"] = {
        "method": "CPM + log2(1+x)",
        "target_sum": float(target_sum),
        "kept_raw_counts_layer": "counts",
        "stored_layer": "log2_cpm"
    }

    # # 3.7 Drive'a KAYDET
    adata.write_h5ad(out_path)

    # # 3.8 Kaydedildi mi? (garanti kontrol)
    saved_ok = os.path.exists(out_path)

    # # 3.9 Kısa özet yazdır
    print("\nSaved file:", out_path)
    print("Saved OK? :", saved_ok)
    print("Shape    :", adata.shape)
    print("Layers   :", list(adata.layers.keys()))

# ----------------------------------------
# 4) Tüm QC2 dosyaları için çalıştır
# ----------------------------------------
for in_path in qc2_paths:
    # # output adı: aynı isim + __log2cpm
    base = os.path.basename(in_path).replace(".h5ad", "")
    out_path = f"{OUT_DIR}/{base}__log2cpm.h5ad"
    normalize_to_log2cpm(in_path, out_path)

print("\nDONE: QC2 normalization completed and saved to Drive.")
