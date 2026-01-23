import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.image as mpimg

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

counts_root = base / "extracted" / "filtered_count_matrices" / "filtered_count_matrices"
spatial_root = base / "extracted" / "spatial" / "spatial"

for sample in samples:
    # --- paths ---
    mtx_dir = counts_root / f"{sample}_filtered_count_matrix"
    sp_dir  = spatial_root / f"{sample}_spatial"
    img_path = sp_dir / "tissue_hires_image.png"
    pos_path = sp_dir / "tissue_positions_list.csv"

    out_dir = base / "outputs" / "STEP06_SPATIAL_QC_HEATMAPS" / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- load counts (fast, and enough for QC) ---
    adata = sc.read_mtx(mtx_dir / "matrix.mtx").T
    barcodes = pd.read_csv(mtx_dir / "barcodes.tsv", header=None, sep="\t")[0].astype(str).values
    adata.obs_names = barcodes

    # spot-level qc
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    qc = adata.obs[["total_counts", "n_genes_by_counts"]].copy()
    qc["barcode"] = qc.index

    # --- load spatial positions + image ---
    pos = pd.read_csv(pos_path, header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    pos_in = pos[pos["in_tissue"] == 1].copy()

    # merge: keep only in-tissue spots that exist in qc
    merged = pos_in.merge(qc, on="barcode", how="inner")

    img = mpimg.imread(img_path)

    # --- heatmap 1: total_counts ---
    plt.figure(figsize=(8,8))
    plt.imshow(img)
    plt.scatter(merged["pxl_col"], merged["pxl_row"], c=merged["total_counts"], s=6, alpha=0.8)
    plt.axis("off")
    plt.title(f"{sample} — Spatial QC: total_counts (in-tissue)")
    plt.colorbar(label="total_counts (UMI/spot)", fraction=0.046, pad=0.04)
    out1 = out_dir / f"{sample}__STEP06__SPATIAL_QC__total_counts.png"
    plt.savefig(out1, dpi=200, bbox_inches="tight")
    plt.close()

    # --- heatmap 2: n_genes_by_counts ---
    plt.figure(figsize=(8,8))
    plt.imshow(img)
    plt.scatter(merged["pxl_col"], merged["pxl_row"], c=merged["n_genes_by_counts"], s=6, alpha=0.8)
    plt.axis("off")
    plt.title(f"{sample} — Spatial QC: n_genes_by_counts (in-tissue)")
    plt.colorbar(label="n_genes_by_counts (genes/spot)", fraction=0.046, pad=0.04)
    out2 = out_dir / f"{sample}__STEP06__SPATIAL_QC__n_genes.png"
    plt.savefig(out2, dpi=200, bbox_inches="tight")
    plt.close()

    print("✅", sample, "| merged spots:", merged.shape[0], "| saved:", out1.name, "and", out2.name)
