import pandas as pd
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

meta_root = base / "extracted" / "metadata" / "metadata"
spatial_root = base / "extracted" / "spatial" / "spatial"

out_dir = base / "outputs" / "STEP06_FILTER_LISTS"
out_dir.mkdir(parents=True, exist_ok=True)

rows = []

for sample in samples:
    meta = pd.read_csv(meta_root / f"{sample}_metadata.csv").rename(columns={"Unnamed: 0":"barcode"})
    meta["barcode"] = meta["barcode"].astype(str)
    meta["is_artefact"] = meta["Classification"].astype(str).str.strip().str.lower().eq("artefact")

    pos = pd.read_csv(spatial_root / f"{sample}_spatial" / "tissue_positions_list.csv", header=None)
    pos.columns = ["barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"]
    pos["barcode"] = pos["barcode"].astype(str)

    pos_in = pos[pos["in_tissue"] == 1][["barcode"]].copy()

    merged = pos_in.merge(meta[["barcode","is_artefact"]], on="barcode", how="left")
    merged["is_artefact"] = merged["is_artefact"].fillna(False).astype(bool)

    keep_main = merged[~merged["is_artefact"]][["barcode"]].copy()
    drop_art  = merged[merged["is_artefact"]][["barcode"]].copy()
    keep_all  = merged[["barcode"]].copy()

    # save lists
    p_keep_main = out_dir / f"{sample}__KEEP_in_tissue_no_artefact.csv"
    p_drop_art  = out_dir / f"{sample}__DROP_artefact.csv"
    p_keep_all  = out_dir / f"{sample}__KEEP_in_tissue_ALL.csv"

    keep_main.to_csv(p_keep_main, index=False)
    drop_art.to_csv(p_drop_art, index=False)
    keep_all.to_csv(p_keep_all, index=False)

    rows.append({
        "sample": sample,
        "in_tissue": int(len(keep_all)),
        "artefact": int(len(drop_art)),
        "keep_main": int(len(keep_main)),
        "paths_saved": "yes"
    })

df = pd.DataFrame(rows).sort_values("sample")
summary_path = out_dir / "STEP06__KEEP_DROP_summary.csv"
df.to_csv(summary_path, index=False)

print("✅ Saved folder:", out_dir)
print("✅ Saved summary:", summary_path)
display(df)

