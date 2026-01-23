from pathlib import Path

base = Path("/content/drive/MyDrive/TÜBİTAK")
root = base / "extracted" / "filtered_count_matrices" / "filtered_count_matrices"
samples = ["1142243F","1160920F","CID4290","CID4465","CID44971","CID4535"]

targets = ["matrix.mtx", "barcodes.tsv", "features.tsv"]

def rename_if_needed(sample):
    sdir = root / f"{sample}_filtered_count_matrix"
    changed = []

    # if plain files already exist, skip
    if all((sdir / t).exists() for t in targets):
        return changed

    # otherwise, rename *.gz -> plain if present
    mapping = {
        "matrix.mtx.gz": "matrix.mtx",
        "barcodes.tsv.gz": "barcodes.tsv",
        "features.tsv.gz": "features.tsv",
    }
    for src, dst in mapping.items():
        src_p = sdir / src
        dst_p = sdir / dst
        if src_p.exists() and (not dst_p.exists()):
            src_p.rename(dst_p)
            changed.append(f"{src} -> {dst}")

    return changed

print("=== Fixing fake .gz filenames (batch) ===")
for sample in samples:
    changes = rename_if_needed(sample)
    sdir = root / f"{sample}_filtered_count_matrix"
    ok = all((sdir / t).exists() for t in targets)
    print(f"\n{sample}:")
    if changes:
        for c in changes: print("  ✅", c)
    else:
        print("  ℹ️ no rename needed (or nothing to rename)")
    print("  ✅ files present:", ok, "|", [t for t in targets if (sdir/t).exists()])
