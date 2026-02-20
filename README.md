Spatial Transcriptomics Analysis of Breast Cancer
TÜBİTAK 2209-A Research Project
Project Overview

This project investigates cancer stem cell (CSC) potential in breast cancer using spatial transcriptomics data.

The main objective is to determine whether the invasive front of the tumor shows different stemness potential compared to the tumor core.

To estimate stemness, I used CytoTRACE2, which predicts differentiation potential based on gene expression diversity.

Research Hypothesis

Cancer progression and metastasis are often driven by stem-like tumor cells.

The hypothesis of this project is:

The invasive front of breast tumors has higher cancer stem cell potential than the tumor core.

To test this hypothesis, I built a systematic spatial analysis pipeline including quality control, deconvolution, cancer-rich filtering, and potency scoring.

Dataset Description
Spatial Transcriptomics Data

Platform: 10x Genomics Visium

Total patients: 6

2 ER+

4 TNBC

Each patient analyzed independently

scRNA-seq Reference Dataset

GEO accession: GSE176078

Used as reference for cell type deconvolution

9 major cell types (based on celltype_major annotation)

Analysis Workflow

The analysis was performed step by step to avoid technical bias and ensure biological relevance.

1. File Validation and Data Integrity

All downloaded files were manually checked.

Metadata, spatial matrices, and filtered matrices were verified.

Some files had incorrect .gz extensions and were fixed.

All samples were successfully loaded after correction.

2. Quality Control (QC)

For each patient, the following metrics were calculated:

total_counts

n_genes_by_counts

mitochondrial gene percentage

Actions performed:

Histogram visualization per patient

Cross-patient QC comparison

1–99 percentile filtering

Removal of high mitochondrial spots

This ensured removal of low-quality and extreme outlier spots.

3. Spatial Validation

Spot coordinates were overlaid on H&E images.

In-tissue spots were visually confirmed.

Artefact regions were identified using the Classification column.

Artefact spots were:

Visualized spatially

Statistically compared to normal spots

Removed from the main dataset

Two filtered sets were created:

Main set: in-tissue + artefact removed

Control set: in-tissue only

Cleaned .h5ad files were saved for each patient.

4. Normalization

For each patient:

Raw counts were preserved.

CPM normalization was applied.

Log2 transformation was generated for visualization.

Raw counts were kept for count-based downstream analysis.

5. Stereoscope Deconvolution

Because Visium spots contain mixed cell populations, direct stemness scoring would be confounded by cell composition.

To solve this:

The scRNA reference dataset was reconstructed into AnnData format.

Matrix orientation errors were corrected.

Shared genes between spatial and scRNA datasets were identified.

Stereoscope (scvi-tools) was used for deconvolution.

Each spatial spot received proportions for 9 major cell types.

Quality control confirmed:

No missing values

Proportions summed to 1 per spot

6. Cancer-Rich Spot Selection

To avoid stromal bias in CytoTRACE2 analysis:

Candidate set: Cancer epithelial proportion > 10%

Strict set: Top 10% within each patient

An adaptive threshold was used because fixed high thresholds removed too many spots due to the mixed-cell nature of Visium.

7. CytoTRACE2 Analysis

CytoTRACE2 was applied to:

CPM-normalized (non-log) expression data

Human species setting

Technical issues encountered were resolved by:

Correct package installation

Ensuring non-log input

Fixing environment inconsistencies

CytoTRACE2 scores were:

Merged back into spatial objects

Visualized on tissue coordinates

Spatial patterns showed non-random clustering.

8. Invasive Front Definition

Instead of using geometric distance, pathology annotations were used.

Definition:

Front: invasive spot adjacent to stroma, adipose, DCIS, or normal tissue

Core: invasive spot surrounded only by invasive spots

This definition is biologically more meaningful than simple center–edge distance.

9. Statistical Analysis

For each patient:

Mann–Whitney U test was performed

Front vs Core CytoTRACE2 scores were compared

FDR correction was applied

Additionally:

All patients were pooled

A global comparison was performed

Results

Patient-level analysis showed heterogeneity.

Some patients had significant front vs core differences, while others did not.

However, pooled analysis across all patients showed:

Invasive front had significantly higher CytoTRACE2 scores

p = 0.0032

This supports the hypothesis that invasive regions may contain more stem-like tumor cells.

Interpretation

The results suggest:

CSC potential is spatially structured

The invasive front may represent a stem-like niche

Breast cancer shows strong inter-patient heterogeneity

Differences between patients may reflect:

Molecular subtype differences

Microenvironment variation

Tumor architecture diversity

Limitations

Visium spots contain mixed cells (not single-cell resolution)

Pathology annotation resolution may affect front definition

Sample size is limited (n = 6)

Future Directions

Planned extensions include:

EMT and proliferation signature analysis

Subtype comparison (ER+ vs TNBC)

Immune microenvironment correlation

Integration with additional spatial features

Tools and Libraries Used

Python

Scanpy

AnnData

scvi-tools (Stereoscope)

CytoTRACE2

Google Colab environment

Conclusion

This project establishes a reproducible spatial transcriptomics pipeline combining:

Quality control

Deconvolution

Cancer-rich filtering

Stemness scoring

Invasive front analysis

The findings support a biologically meaningful increase in stem-like potential at the invasive tumor front.
