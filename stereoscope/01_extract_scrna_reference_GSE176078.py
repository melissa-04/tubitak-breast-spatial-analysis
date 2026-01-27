# ============================================
# STEP09B — Extract scRNA reference (GSE176078)
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, tarfile

ROOT = "/content/drive/MyDrive/TÜBİTAK"
SCRNA_TGZ = f"{ROOT}/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

OUT_DIR = f"{ROOT}/reference/GSE176078/scrna"
os.makedirs(OUT_DIR, exist_ok=True)

print("scRNA tgz exists?:", os.path.exists(SCRNA_TGZ))
print("Extract to:", OUT_DIR)

# 1) İçeriği listele (ilk 60)
with tarfile.open(SCRNA_TGZ, "r:gz") as tar:
    names = tar.getnames()
    print("\n--- scRNA archive content (first 60) ---")
    for n in names[:60]:
        print(" -", n)
    print("Total files:", len(names))

# 2) Çıkar (tek seferlik marker ile)
marker = f"{OUT_DIR}/__EXTRACTED_OK__"
if not os.path.exists(marker):
    print("\nExtracting now...")
    with tarfile.open(SCRNA_TGZ, "r:gz") as tar:
        tar.extractall(path=OUT_DIR)
    open(marker, "w").write("ok")
    print("Extraction done ✅")
else:
    print("\nAlready extracted ✅")
