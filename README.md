# TCGA Pan-Cancer Gene Expression Clustering

A lightweight Python pipeline performing unsupervised dimensionality reduction (PCA) on transcriptomic gene expression profiles across multiple cancer types.

## Overview
This repository demonstrates how high-dimensional gene expression data can be reduced down to two principal components (PC1 and PC2) to reveal distinct clustering across four carcinoma subtypes:
- **BRCA** (Breast)
- **LUAD** (Lung)
- **PRAD** (Prostate)
- **KIRC** (Kidney)

## Running the Script
```bash
python3 tcga_analysis.py

git add tcga_analysis.py tcga_pca_plot.png README.md
git commit -m "Initial commit: TCGA PCA clustering analysis"
