import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.image as mpimg

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

meta_root = base / "extracted" / "metadata" / "metadata"
spatial_root = base / "extracted" / "spatial" / "spatial"

for sample in samples:
    meta_path = meta_root / f"{sample}_metadata.csv"
    sp_dir = spatial_root / f"{sample}_spatial"
    pos_path = sp_dir / "tissue_positions_list.csv"
    img_path = sp_dir / "tissue_hires_image.png"

    out_dir = base / "outputs" / "STEP06_ARTEFACT_MAP" / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / f"{sample}__STEP06__ARTEFACT_map.png"

    # read metadata
    meta = pd.read_csv(meta_path).rename(columns={"Unnamed: 0": "barcode"})
    meta["barcode"] = meta["barcode"].astype(str)
    meta["is_artefact"] = meta["Classification"].astype(str).str.strip().str.lower().eq("artefact")

    # read positions
    pos = pd.read_csv(pos_path, header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    pos_in = pos[pos["in_tissue"] == 1].copy()
    pos_in["barcode"] = pos_in["barcode"].astype(str)

    # merge
    merged = pos_in.merge(meta[["barcode","is_artefact"]], on="barcode", how="left")
    merged["is_artefact"] = merged["is_artefact"].fillna(False)

    n_total = merged.shape[0]
    n_art = int(merged["is_artefact"].sum())

    # plot
    img = mpimg.imread(img_path)
    plt.figure(figsize=(8,8))
    plt.imshow(img)

    normal = merged[~merged["is_artefact"]]
    art = merged[merged["is_artefact"]]

    plt.scatter(normal["pxl_col"], normal["pxl_row"], s=5, alpha=0.55)
    plt.scatter(art["pxl_col"], art["pxl_row"], s=12, alpha=0.95, c="red")

    plt.axis("off")
    plt.title(f"{sample} — Artefact map (red={n_art}, in-tissue={n_total})")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()

    print(f"✅ {sample} | in-tissue={n_total} | artefact={n_art} | saved: {out_png}")
