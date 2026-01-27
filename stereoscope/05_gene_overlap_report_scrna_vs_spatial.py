# ============================================
# STEP09G — Gene overlap report (scRNA vs spatial)
# Saves CSV report to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # scRNA reference (with cell_type)
SCRNA_PATH = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE/GSE176078_scRNA_reference_counts__label_celltype_major.h5ad"

# # spatial counts prepared for stereoscope
SPATIAL_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SPATIAL_QC2_COUNTS"

# # report output
REPORT_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/REPORTS"
os.makedirs(REPORT_DIR, exist_ok=True)
OUT_CSV = f"{REPORT_DIR}/gene_overlap_report.csv"

# # load scRNA genes
scrna = sc.read_h5ad(SCRNA_PATH, backed="r")
scrna_genes = set(scrna.var_names)

# # list spatial files
spatial_files = sorted([f for f in os.listdir(SPATIAL_DIR) if f.endswith("__counts_for_stereoscope.h5ad")])

print("scRNA genes:", len(scrna_genes))
print("Spatial files:", len(spatial_files))

rows = []
for f in spatial_files:
    sp = sc.read_h5ad(f"{SPATIAL_DIR}/{f}", backed="r")
    sp_genes = set(sp.var_names)

    overlap = scrna_genes.intersection(sp_genes)

    rows.append({
        "spatial_file": f,
        "spatial_genes": len(sp_genes),
        "scrna_genes": len(scrna_genes),
        "overlap_genes": len(overlap)
    })

    print(f"\n{f}")
    print(" - spatial genes:", len(sp_genes))
    print(" - overlap genes:", len(overlap))

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)

print("\nSAVED TO DRIVE ✅")
print("Path:", OUT_CSV)
print("Saved OK?:", os.path.exists(OUT_CSV))
print("\nPreview:")
print(df)
