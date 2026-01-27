# ============================================
# STEP 01 — Setup environment and list h5ad inputs
# ============================================

# 1) Google Drive'ı bağla
from google.colab import drive
drive.mount('/content/drive')

# 2) Proje ana klasörünü tanımla
# Burayı değiştirme: Drive ekranında gösterdiğin yapıya göre ayarlı
PROJECT_ROOT = "/content/drive/MyDrive/TÜBİTAK"

# 3) QC'den geçmiş h5ad dosyalarının olduğu klasör
# (STEP07B_H5AD_QC_FILTERED içindekiler bizim normalizasyon girdimiz)
H5AD_INPUT_DIR = f"{PROJECT_ROOT}/outputs/STEP07B_H5AD_QC_FILTERED"

# 4) Klasör var mı kontrol et
import os
print("Klasör mevcut mu?:", os.path.exists(H5AD_INPUT_DIR))
print("Klasör yolu:", H5AD_INPUT_DIR)

# 5) İçindeki .h5ad dosyalarını listele
h5ad_files = [f for f in os.listdir(H5AD_INPUT_DIR) if f.endswith(".h5ad")]

print("\nBulunan h5ad dosyaları:")
for f in sorted(h5ad_files):
    print(" -", f)

print(f"\nToplam h5ad sayısı: {len(h5ad_files)}")
