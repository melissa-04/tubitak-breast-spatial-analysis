import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
summary_path = base / "outputs" / "STEP06_QC" / "STEP06__QC_summary_ALL_samples.csv"

df = pd.read_csv(summary_path).sort_values("sample")

out_dir = base / "outputs" / "STEP06_QC"
out_png1 = out_dir / "STEP06__QC_ALL__median_total_counts_bar.png"
out_png2 = out_dir / "STEP06__QC_ALL__median_n_genes_bar.png"

# Plot 1: median_total_counts
plt.figure(figsize=(8,4))
plt.bar(df["sample"], df["median_total_counts"])
plt.xlabel("Sample")
plt.ylabel("Median total_counts (UMI/spot)")
plt.title("QC summary — Median total counts per spot (all samples)")
plt.xticks(rotation=30, ha="right")
plt.savefig(out_png1, dpi=200, bbox_inches="tight")
plt.close()

# Plot 2: median_n_genes_by_counts
plt.figure(figsize=(8,4))
plt.bar(df["sample"], df["median_n_genes_by_counts"])
plt.xlabel("Sample")
plt.ylabel("Median n_genes_by_counts (genes/spot)")
plt.title("QC summary — Median genes detected per spot (all samples)")
plt.xticks(rotation=30, ha="right")
plt.savefig(out_png2, dpi=200, bbox_inches="tight")
plt.close()

print("✅ Read:", summary_path)
print("✅ Saved:", out_png1)
print("✅ Saved:", out_png2)
display(df)
