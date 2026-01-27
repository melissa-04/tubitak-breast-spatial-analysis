# ============================================
# STEP10E — Run SpatialStereoscope for ALL samples
# Saves per-sample CSV + h5ad + a run summary to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os, re
import pandas as pd
import scanpy as sc
from scvi.external import RNAStereoscope, SpatialStereoscope

ROOT = "/content/drive/MyDrive/TÜBİTAK"

# # Inputs
MODEL_DIR   = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/MODELS/rna_model"
SCRNA_TRAIN = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/INPUTS/scrna_train_downsampled.h5ad"
INPUT_DIR   = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/INPUTS"

# # Outputs
RESULTS_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/RESULTS"
REPORT_DIR  = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/REPORTS"
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)

# # Helper: safe column name
def safe_name(x):
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    return x.strip("_")

# # 1) Load scRNA train + load RNA model
adata_sc = sc.read_h5ad(SCRNA_TRAIN)
adata_sc.layers["counts"] = adata_sc.X.copy()
RNAStereoscope.setup_anndata(adata_sc, layer="counts", labels_key="cell_type")
rna_model = RNAStereoscope.load(MODEL_DIR, adata=adata_sc)
print("RNA model loaded ✅")

# # 2) Find spatial intersect inputs
spatial_files = sorted([f for f in os.listdir(INPUT_DIR) if f.startswith("spatial_") and f.endswith("_counts_intersect.h5ad")])
print("Spatial inputs found:", len(spatial_files))
for f in spatial_files:
    print(" -", f)

summary_rows = []

# # 3) Loop over samples
for f in spatial_files:
    sample = f.replace("spatial_", "").replace("_counts_intersect.h5ad", "")
    in_path = f"{INPUT_DIR}/{f}"

    out_dir = f"{RESULTS_DIR}/{sample}"
    os.makedirs(out_dir, exist_ok=True)

    out_csv  = f"{out_dir}/stereoscope_proportions_{sample}.csv"
    out_h5ad = f"{out_dir}/spatial_{sample}_stereoscope_result.h5ad"

    print(f"\n=== Running sample: {sample} ===")
    print("Input:", in_path)

    adata_sp = sc.read_h5ad(in_path)
    adata_sp.layers["counts"] = adata_sp.X.copy()
    SpatialStereoscope.setup_anndata(adata_sp, layer="counts")

    sp_model = SpatialStereoscope.from_rna_model(adata_sp, rna_model)
    sp_model.train(max_epochs=200)

    props = sp_model.get_proportions()
    if not isinstance(props, pd.DataFrame):
        props = pd.DataFrame(props)

    # # Save CSV
    props.to_csv(out_csv, index=True)

    # # Save into h5ad (obs + obsm)
    adata_sp.obsm["stereoscope_proportions"] = props.values
    adata_sp.uns["stereoscope_cell_types"] = list(props.columns)
    for col in props.columns:
        adata_sp.obs[f"stereo_{safe_name(col)}"] = props[col].values

    adata_sp.write_h5ad(out_h5ad)

    # # Summary
    summary_rows.append({
        "sample": sample,
        "n_spots": adata_sp.n_obs,
        "n_genes": adata_sp.n_vars,
        "n_cell_types": props.shape[1],
        "csv_path": out_csv,
        "h5ad_path": out_h5ad
    })

    print("Saved OK?:", os.path.exists(out_csv), os.path.exists(out_h5ad))

# # 4) Save run summary
summary_df = pd.DataFrame(summary_rows)
summary_path = f"{REPORT_DIR}/stereoscope_run_summary.csv"
summary_df.to_csv(summary_path, index=False)

print("\nSAVED TO DRIVE ✅ (run summary)")
print("Path:", summary_path)
print("Saved OK?:", os.path.exists(summary_path))
print(summary_df)
print("\nDONE ✅")
