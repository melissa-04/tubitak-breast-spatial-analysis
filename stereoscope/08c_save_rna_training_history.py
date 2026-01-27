# ============================================
# STEP10C-FIX — Save RNA model training history safely
# Saves CSV to Google Drive
# ============================================

from google.colab import drive
drive.mount('/content/drive')

import os
import pandas as pd

ROOT = "/content/drive/MyDrive/TÜBİTAK"
REPORT_DIR = f"{ROOT}/outputs/STEP10_STEREOSCOPE_RUN/REPORTS"
os.makedirs(REPORT_DIR, exist_ok=True)

# # sc_model.history farklı formatlarda gelebilir, hepsini güvenli yakala
hist_obj = sc_model.history

if hasattr(hist_obj, "to_dataframe"):
    # # bazı sürümlerde History objesi böyle gelir
    hist_df = hist_obj.to_dataframe()
elif isinstance(hist_obj, pd.DataFrame):
    hist_df = hist_obj
elif isinstance(hist_obj, dict):
    # # dict values list ise DataFrame(hist_obj) olur; scalar ise tek satır yap
    try:
        hist_df = pd.DataFrame(hist_obj)
    except Exception:
        hist_df = pd.DataFrame([hist_obj])
else:
    # # en kötü ihtimal: tip bilgisini kaydet
    hist_df = pd.DataFrame({"history_type": [str(type(hist_obj))]})

out_path = f"{REPORT_DIR}/rna_model_training_history.csv"
hist_df.to_csv(out_path, index=False)

print("SAVED TO DRIVE ✅ (training history)")
print("Path:", out_path)
print("Saved OK?:", os.path.exists(out_path))
print("History shape:", hist_df.shape)
print(hist_df.head())
