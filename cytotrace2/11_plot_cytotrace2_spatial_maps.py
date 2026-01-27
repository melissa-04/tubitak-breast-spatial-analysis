# ============================================
# STEP14 — Plot CytoTRACE2 spatial maps (strict cancer-rich spots)
# Saves PNG figures to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import scanpy as sc
import matplotlib.pyplot as plt

ROOT = "/content/drive/MyDrive/TÜBİTAK"
IN_BASE = f"{ROOT}/outputs/STEP13_CYTOTRACE2/MERGED"
OUT_DIR = f"{ROOT}/outputs/STEP14_FIGURES/CYTOTRACE2_SPATIAL"
os.makedirs(OUT_DIR, exist_ok=True)

samples = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
KEY = "cytotrace2_potency_score"

for s in samples:
    in_h5ad = f"{IN_BASE}/{s}/spatial_{s}__strict__cytotrace2_scored.h5ad"
    adata = sc.read_h5ad(in_h5ad)

    out_png = f"{OUT_DIR}/cytotrace2_spatial_{s}.png"

    # Try Scanpy's spatial plot (requires spatial info)
    ok_spatial = ("spatial" in adata.obsm_keys()) and ("spatial" in adata.uns_keys() or True)

    plt.figure()
    try:
        sc.pl.spatial(adata, color=KEY, show=False)
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"{s}: saved spatial plot ✅", os.path.exists(out_png))
    except Exception as e:
        # Fallback: simple scatter using coordinates
        if "spatial" in adata.obsm_keys():
            coords = adata.obsm["spatial"]
            plt.scatter(coords[:,0], coords[:,1], c=adata.obs[KEY].values, s=8)
            plt.gca().invert_yaxis()
            plt.title(f"{s} — CytoTRACE2 potency (strict)")
            plt.colorbar()
            plt.savefig(out_png, dpi=200, bbox_inches="tight")
            plt.close()
            print(f"{s}: fallback plot ✅", os.path.exists(out_png))
        else:
            print(f"{s}: ERROR no spatial coords; cannot plot. Error:", str(e))

print("\nDONE ✅")
print("Figures saved to:", OUT_DIR)
