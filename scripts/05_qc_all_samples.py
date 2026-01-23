import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
root = base / "extracted" / "filtered_count_matrices" / "filtered_count_matrices"

samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]
all_rows = []

for sample in samples:
    mtx_dir = root / f"{sample}_filtered_count_matrix"
    out_dir = base / "outputs" / "STEP06_QC" / sample
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- CLEAN old outputs (only this step's files) ---
    for old in out_dir.glob(f"{sample}__STEP06__*"):
        old.unlink()

    # --- load ---
    adata = sc.read_mtx(mtx_dir / "matrix.mtx").T

    genes = pd.read_csv(mtx_dir / "features.tsv", header=None, sep="\t")
    barcodes = pd.read_csv(mtx_dir / "barcodes.tsv", header=None, sep="\t")
    adata.var_names = genes.iloc[:,0].astype(str).values
    adata.obs_names = barcodes.iloc[:,0].astype(str).values
    adata.var_names_make_unique()

    # --- QC ---
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    row = {
        "sample": sample,
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "median_total_counts": float(adata.obs["total_counts"].median()),
        "median_n_genes_by_counts": float(adata.obs["n_genes_by_counts"].median()),
    }
    all_rows.append(row)

    # --- save per-sample summary ---
    pd.DataFrame([row]).to_csv(out_dir / f"{sample}__STEP06__QC_summary.csv", index=False)

    # --- plots ---
    plt.figure()
    adata.obs["total_counts"].hist(bins=50)
    plt.xlabel("Total counts (UMI per spot)")
    plt.ylabel("Number of spots")
    plt.title(f"{sample} — Total counts per spot")
    plt.savefig(out_dir / f"{sample}__STEP06__HIST_total_counts.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure()
    adata.obs["n_genes_by_counts"].hist(bins=50)
    plt.xlabel("Number of genes per spot")
    plt.ylabel("Number of spots")
    plt.title(f"{sample} — Genes detected per spot")
    plt.savefig(out_dir / f"{sample}__STEP06__HIST_n_genes.png", dpi=200, bbox_inches="tight")
    plt.close()

    print("✅ QC done:", sample)

# --- combined summary ---
all_df = pd.DataFrame(all_rows).sort_values("sample")
combined_path = base / "outputs" / "STEP06_QC" / "STEP06__QC_summary_ALL_samples.csv"
all_df.to_csv(combined_path, index=False)

print("\n✅ Saved combined summary:", combined_path)
display(all_df)
