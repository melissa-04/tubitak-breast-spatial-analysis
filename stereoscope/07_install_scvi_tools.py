# ============================================
# STEP10B — Install scvi-tools (Stereoscope package)
# ============================================

# # Stereoscope, scvi-tools içinde geliyor
!pip -q install "scvi-tools>=1.2.0"

# # Sürüm ve GPU kontrolü
import scvi, torch
print("scvi-tools:", scvi.__version__)
print("torch:", torch.__version__)
print("GPU available?:", torch.cuda.is_available())
