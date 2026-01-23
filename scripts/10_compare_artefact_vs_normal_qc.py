import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

counts_root = base / "extracted" / "filtered_count_matrices" / "filtered_count_matrices"
meta_root   = base / "extracted" / "metadata" / "metadata"
spatial_root= base / "extracted" / "spatial" / "spatial"

out_dir = base / "outputs" / "STEP06_ARTEFACT_QC_COMPARE"
out_dir.mkdir(parents=True, exist_ok=True)

rows = []

for sample in samples:
    # --- load counts (only need barcodes + qc) ---
    mtx_dir = counts_root / f"{sample}_filtered_count_matrix"
    adata = sc.read_mtx(mtx_dir / "matrix.mtx").T
    barcodes = pd.read_csv(mtx_dir / "barcodes.tsv", header=None, sep="\t")[0].astype(str).values
    adata.obs_names = barcodes
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    qc = adata.obs[["total_counts", "n_genes_by_counts"]].copy()
    qc["barcode"] = qc.index.astype(str)

    # --- load metadata (artefact flag) ---
    meta = pd.read_csv(meta_root / f"{sample}_metadata.csv").rename(columns={"Unnamed: 0":"barcode"})
    meta["barcode"] = meta["barcode"].astype(str)
    meta["is_artefact"] = meta["Classification"].astype(str).str.strip().str.lower().eq("artefact")

    # --- load spatial positions (in-tissue only) ---
    pos = pd.read_csv(spatial_root / f"{sample}_spatial" / "tissue_positions_list.csv", header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    pos = pos[pos["in_tissue"] == 1].copy()
    pos["barcode"] = pos["barcode"].astype(str)

    # --- merge: in-tissue + qc + artefact ---
    merged = pos[["barcode"]].merge(qc, on="barcode", how="inner").merge(meta[["barcode","is_artefact"]], on="barcode", how="left")
    merged["is_artefact"] = merged["is_artefact"].fillna(False).astype(bool)  # FutureWarning fix

    art = merged[merged["is_artefact"]]
    norm = merged[~merged["is_artefact"]]

    row = {
        "sample": sample,
        "n_in_tissue": int(merged.shape[0]),
        "n_artefact": int(art.shape[0]),
        "artefact_rate_%": float(100 * art.shape[0] / max(1, merged.shape[0])),
        "median_total_counts_normal": float(norm["total_counts"].median()) if len(norm) else float("nan"),
        "median_total_counts_artefact": float(art["total_counts"].median()) if len(art) else float("nan"),
        "median_n_genes_normal": float(norm["n_genes_by_counts"].median()) if len(norm) else float("nan"),
        "median_n_genes_artefact": float(art["n_genes_by_counts"].median()) if len(art) else float("nan"),
    }
    rows.append(row)

df = pd.DataFrame(rows).sort_values("sample")
csv_path = out_dir / "STEP06__artefact_vs_normal_qc_summary.csv"
df.to_csv(csv_path, index=False)

print("✅ Saved table:", csv_path)
display(df)

# --- Plot 1: median total_counts (normal vs artefact) ---
plt.figure(figsize=(9,4))
x = range(len(df))
plt.bar([i - 0.2 for i in x], df["median_total_counts_normal"], width=0.4, label="normal")
plt.bar([i + 0.2 for i in x], df["median_total_counts_artefact"], width=0.4, label="artefact")
plt.xticks(list(x), df["sample"], rotation=30, ha="right")
plt.ylabel("Median total_counts (UMI/spot)")
plt.title("Artefact vs Normal — Median total_counts")
plt.legend()
p1 = out_dir / "STEP06__artefact_vs_normal__median_total_counts.png"
plt.savefig(p1, dpi=200, bbox_inches="tight")
plt.close()

# --- Plot 2: median n_genes ---
plt.figure(figsize=(9,4))
plt.bar([i - 0.2 for i in x], df["median_n_genes_normal"], width=0.4, label="normal")
plt.bar([i + 0.2 for i in x], df["median_n_genes_artefact"], width=0.4, label="artefact")
plt.xticks(list(x), df["sample"], rotation=30, ha="right")
plt.ylabel("Median n_genes_by_counts (genes/spot)")
plt.title("Artefact vs Normal — Median n_genes")
plt.legend()
p2 = out_dir / "STEP06__artefact_vs_normal__median_n_genes.png"
plt.savefig(p2, dpi=200, bbox_inches="tight")
plt.close()

print("✅ Saved plots:", p1)
print("✅ Saved plots:", p2)
