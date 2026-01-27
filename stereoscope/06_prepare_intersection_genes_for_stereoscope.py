# ============================================
# STEP10A — Prepare intersected gene set for Stereoscope
# Saves intersected inputs to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import scanpy as sc

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # scRNA reference (counts + cell_type)
SCRNA_PATH = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SCRNA_REFERENCE/GSE176078_scRNA_reference_counts__label_celltype_major.h5ad"

# # spatial counts for stereoscope
SPATIAL_DIR = f"{ROOT}/outputs/STEP09_STEREOSCOPE_PREP/SPATIAL_QC2_COUNTS"
spatial_files = sorted([f for f in os.listdir(SPATIAL_DIR) if f.endswith("__counts_for_stereoscope.h5ad")])

# # output dirs
OUT_BASE = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN"
IN_DIR   = f"{OUT_BASE}/INPUTS"
REP_DIR  = f"{OUT_BASE}/REPORTS"
os.makedirs(IN_DIR, exist_ok=True)
os.makedirs(REP_DIR, exist_ok=True)

print("Spatial files:", len(spatial_files))

# # 1) Load scRNA (backed ok) and start intersection set
scrna = sc.read_h5ad(SCRNA_PATH)
gene_set = set(scrna.var_names)

# # 2) Intersect with each spatial gene list
for f in spatial_files:
    sp = sc.read_h5ad(f"{SPATIAL_DIR}/{f}", backed="r")
    gene_set = gene_set.intersection(set(sp.var_names))

gene_list = sorted(gene_set)
print("Final intersection gene count:", len(gene_list))

# # Save gene list
gene_txt = f"{REP_DIR}/intersection_gene_list.txt"
with open(gene_txt, "w") as w:
    for g in gene_list:
        w.write(g + "\n")

print("\nSAVED gene list ✅")
print("Path:", gene_txt)
print("Saved OK?:", os.path.exists(gene_txt))

# # 3) Subset and save scRNA
scrna_sub = scrna[:, gene_list].copy()
scrna_out = f"{IN_DIR}/scrna_counts_intersect.h5ad"
scrna_sub.write_h5ad(scrna_out)

print("\nSAVED scRNA intersect ✅")
print("Path:", scrna_out)
print("Saved OK?:", os.path.exists(scrna_out))
print("Shape:", scrna_sub.shape)

# # 4) Subset and save each spatial
for f in spatial_files:
    sp = sc.read_h5ad(f"{SPATIAL_DIR}/{f}")
    sp_sub = sp[:, gene_list].copy()

    sample = f.split("__")[0]  # e.g., 1142243F
    sp_out = f"{IN_DIR}/spatial_{sample}_counts_intersect.h5ad"
    sp_sub.write_h5ad(sp_out)

    print("\nSAVED spatial intersect ✅", sample)
    print("Path:", sp_out)
    print("Saved OK?:", os.path.exists(sp_out))
    print("Shape:", sp_sub.shape)

print("\nDONE: Intersected inputs prepared.")
