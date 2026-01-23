import tarfile
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
archive = base / "spatial.tar.gz"
out_dir = base / "extracted" / "spatial"
out_dir.mkdir(parents=True, exist_ok=True)

def safe_extract(tar, path: Path):
    path = path.resolve()
    for m in tar.getmembers():
        mp = (path / m.name).resolve()
        if not str(mp).startswith(str(path)):
            raise Exception("Unsafe path detected in tar!")
    tar.extractall(path)

with tarfile.open(archive, "r:gz") as tar:
    safe_extract(tar, out_dir)

print("✅ Extracted to:", out_dir)

# quick scan: which samples have the key Visium spatial files?
needed_any = [
    "scalefactors_json.json",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
    "tissue_positions_list.csv",
    "tissue_positions.csv"
]

# find candidate sample folders by searching for scalefactors files
scales = list(out_dir.rglob("scalefactors_json.json"))
print("\nFound scalefactors files:", len(scales))
for s in scales[:20]:
    print(" -", s.relative_to(out_dir))

print("\n=== Per-sample check (first 15 hits) ===")
hits = 0
for sf in scales:
    sample_dir = sf.parent
    rel = sample_dir.relative_to(out_dir)
    present = {n: (sample_dir / n).exists() for n in needed_any}
    print("\nSample folder:", rel)
    for k,v in present.items():
        if v: print("  ✅", k)
    # show if tissue positions exists in either name
    if not (present["tissue_positions_list.csv"] or present["tissue_positions.csv"]):
        print("  ❗ No tissue_positions*.csv found in this folder")
    hits += 1
    if hits >= 15:
        break
