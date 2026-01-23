import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.image as mpimg

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

for sample in samples:
    sp_dir = base / "extracted" / "spatial" / "spatial" / f"{sample}_spatial"
    img_path = sp_dir / "tissue_hires_image.png"
    pos_path = sp_dir / "tissue_positions_list.csv"

    out_dir = base / "outputs" / "STEP06_SPATIAL_CHECK" / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{sample}__STEP06__tissue_with_spots.png"

    img = mpimg.imread(img_path)

    pos = pd.read_csv(pos_path, header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    pos_in = pos[pos["in_tissue"] == 1]

    plt.figure(figsize=(8,8))
    plt.imshow(img)
    plt.scatter(pos_in["pxl_col"], pos_in["pxl_row"], s=4, alpha=0.6)
    plt.axis("off")
    plt.title(f"{sample} — tissue_hires + in-tissue spots ({len(pos_in)})")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()

    print("✅ Saved:", out_png)
