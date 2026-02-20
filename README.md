Spatial Transcriptomics Analysis of Breast Cancer

TÜBİTAK 2209-A Project

1. Project Overview

This project investigates cancer stem cell (CSC) potential in breast cancer using spatial transcriptomics data.

The main goal is to understand whether the invasive front of the tumor shows different stemness potential compared to the tumor core.

To measure stemness, I used CytoTRACE2, which estimates differentiation potential based on gene expression profiles.

2. Research Hypothesis

Tumor progression and metastasis are often driven by stem-like cancer cells.

The hypothesis of this project is:

The invasive front of breast tumors has higher cancer stem cell potential compared to the tumor core.

To test this, spatial transcriptomics data was analyzed step by step, including quality control, deconvolution, cancer-rich filtering, and potency scoring.

3. Dataset Description
Spatial Data

Platform: 10x Genomics Visium

6 breast cancer patients

2 ER+

4 TNBC

Each patient was analyzed separately to preserve biological differences.

scRNA-seq Reference

Dataset: GSE176078

Used as reference for cell type deconvolution

9 major cell types used (based on celltype_major annotation)

4. Analysis Workflow

The analysis was performed systematically to avoid technical bias and ensure biological interpretation.

Step 1 — File Validation and Data Integrity

All downloaded files were checked manually:

Metadata

Spatial count matrices

Filtered matrices

Some files had incorrect .gz extensions. These were corrected to ensure proper loading.

All samples were successfully loaded after validation.

Step 2 — Quality Control (QC)

For each patient, I calculated:

total_counts

n_genes_by_counts

mitochondrial gene percentage

I created:

Histograms per patient

Cross-patient QC comparisons

Percentile-based filtering (1–99%)

Removal of high mitochondrial spots

This ensured that low-quality and extreme outlier spots were removed before further analysis.

Step 3 — Spatial Validation

To confirm data integrity:

Spot coordinates were overlaid on H&E images

In-tissue spots were visually verified

Artefact regions were identified using the "Classification" column

Artefact spots were:

Visualized spatially

Compared statistically with normal spots

Removed from the main dataset

Two filtered sets were created:

Main set: in-tissue + artefact removed

Control set: in-tissue only

Cleaned .h5ad files were saved for each patient.

Step 4 — Normalization

For each sample:

Raw counts were preserved

CPM normalization was applied

Log2 transformation was generated for visualization

This step prepared the data for downstream analysis while preserving raw counts for count-based methods.

Step 5 — Stereoscope Deconvolution

Spatial spots contain mixed cell populations. Therefore, direct stemness analysis would be confounded by cell composition.

To address this:

The scRNA reference dataset was reconstructed into AnnData format.

Matrix orientation errors were corrected.

Shared genes between spatial and scRNA datasets were identified.

Stereoscope (scvi-tools) was used for deconvolution.

Each spatial spot received proportions for 9 major cell types.

Quality control confirmed:

No missing values

Proportions summed to 1 per spot

Step 6 — Cancer-Rich Spot Selection

CytoTRACE2 should not be applied to mixed stromal regions.

Therefore, I filtered spots based on cancer epithelial proportion:

Candidate set: Cancer epithelial > 10%

Strict set: Top 10% within each patient

This adaptive threshold was chosen because fixed high thresholds removed too many spots due to Visium’s mixed-cell nature.

Step 7 — CytoTRACE2 Analysis

CytoTRACE2 was applied to:

CPM-normalized (non-log) expression matrices

Human species setting

Initial technical errors were resolved by:

Correct package installation

Ensuring non-log input

Fixing environment inconsistencies

CytoTRACE2 scores were:

Merged back into spatial data

Visualized spatially on tissue coordinates

The scores showed non-random spatial clustering.

Step 8 — Invasive Front Definition

Instead of using geometric distance, I used pathology annotations.

Definition:

Front: invasive spot adjacent to stroma, adipose, DCIS, or normal tissue

Core: invasive spot surrounded only by invasive spots

This method is biologically more meaningful than center–edge distance.

Step 9 — Statistical Analysis

For each patient:

Mann–Whitney U test was performed

Front vs Core CytoTRACE2 scores compared

FDR correction applied

Additionally:

All patients were pooled

Global comparison was performed

5. Results

Patient-level results showed heterogeneity.

Some patients had significant differences between front and core, others did not.

However, pooled analysis across all patients showed:

Invasive front had significantly higher CytoTRACE2 scores

p = 0.0032

This supports the hypothesis that invasive regions may harbor more stem-like tumor cells.

6. Interpretation

The results suggest:

CSC potential is spatially structured

Invasive front may represent a biologically active stem-like niche

Inter-patient heterogeneity is strong in breast cancer

The variability between patients may reflect:

Molecular subtype differences

Microenvironment variation

Tumor architecture differences

7. Limitations

Visium spots contain mixed cells (not single-cell resolution)

Pathology annotation resolution may influence front definition

Sample size is limited (6 patients)
