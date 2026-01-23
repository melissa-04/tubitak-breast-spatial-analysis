import scanpy as sc
import pandas as pd
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

counts_root = base / "extracted" / "filtered_count_matrices" / "filtered_count_matrices"
lists_root  = base / "outputs" / "STEP06_FILTER_LISTS"
out_root    = base / "outputs" / "STEP07_H5AD_MAIN"
out_root.mkdir(parents=True, exist_ok=True)

for sample in samples:
    mtx_dir = counts_root / f"{sample}_filtered_count_matrix"
    keep_path = lists_root / f"{sample}__KEEP_in_tissue_no_artefact.csv"

    out_dir = out_root / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = out_dir / f"{sample}__MAIN_keep_in_tissue_no_artefact.h5ad"

    # load full counts
    adata = sc.read_mtx(mtx_dir / "matrix.mtx").T
    genes = pd.read_csv(mtx_dir / "features.tsv", header=None, sep="\t")
    barcodes = pd.read_csv(mtx_dir / "barcodes.tsv", header=None, sep="\t")[0].astype(str).values

    adata.var_names = genes.iloc[:,0].astype(str).values
    adata.obs_names = barcodes
    adata.var_names_make_unique()

    # load keep list
    keep_raw = pd.read_csv(keep_path)["barcode"].astype(str).tolist()

    # FIX: intersect with existing barcodes in count matrix
    obs_set = set(adata.obs_names)
    keep = [b for b in keep_raw if b in obs_set]
    missing = len(keep_raw) - len(keep)

    if missing > 0:
        print(f"⚠️ {sample}: {missing} barcodes in KEEP list are missing in count matrix -> dropped")

    # subset safely
    adata = adata[keep, :].copy()

    # qc metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # save
    adata.write_h5ad(out_h5ad)
    print(f"✅ Saved {sample}: spots={adata.n_obs}, genes={adata.n_vars} -> {out_h5ad}")

print("\n✅ Done.")
