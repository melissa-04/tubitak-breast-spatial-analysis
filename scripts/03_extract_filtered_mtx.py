import tarfile
from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
archive = base / "filtered_count_matrices (1).tar.gz"
out_dir = base / "extracted" / "filtered_count_matrices"
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

# list sample folders and show whether expected files exist
expected = ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]

# find folders that contain matrix.mtx.gz
matrices = list(out_dir.rglob("matrix.mtx.gz"))
print("\nFound matrix.mtx.gz files:", len(matrices))
for m in matrices[:10]:
    print(" -", m.relative_to(out_dir))

print("\n=== Check first 6 sample folders ===")
seen = 0
for m in matrices:
    sample_dir = m.parent
    rel = sample_dir.relative_to(out_dir)
    print("\nSample folder:", rel)
    for e in expected:
        print("  ", "✅" if (sample_dir/e).exists() else "❌", e)
    seen += 1
    if seen >= 6:
        break
