from google.colab import drive
drive.mount('/content/drive')

import os
from pathlib import Path

root = Path("/content/drive/MyDrive/TÜBİTAK")

if not root.exists():
    raise FileNotFoundError(f"Folder not found: {root}")

print("📁 Folder:", root)
print("\nFiles:\n" + "-"*60)

items = sorted(root.iterdir(), key=lambda p: p.name.lower())
for p in items:
    if p.is_file():
        size_mb = p.stat().st_size / (1024**2)
        print(f"{p.name:35s}  {size_mb:10.2f} MB")
    else:
        print(f"[DIR] {p.name}")

print("-"*60)
