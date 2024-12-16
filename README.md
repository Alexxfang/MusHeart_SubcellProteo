# MusHeart_SubcellProteo

# Script for Generating Outputs Corresponding to Supplementary Tables 12 to 17

## Overview
This script processes proteomics data to generate outputs that correspond to **Supplementary Tables 12 to 17**. These tables summarize the results of machine learning (ML)-based subcellular predictions and protein-protein interaction (PPI) analyses using STRINGdb and GO-based subcellular annotations.

## Prerequisites
Ensure the following files are in the working directory:
1. **Input Data Files**:
   - `norm_transform_lopit_18092023.csv`: 0-1 scaled protein intensity ratio.
   - `Maker_list_25112023_mpv2.csv`: Subcellular Protein markers.
   - `metadata_29082023.csv`: Sample information.
   - `f1.csv`: Expression matrix 0-1 scaled protein intensity ratio
   - `f2.csv`: Feature metadata including marker labels.
   - `f3.csv`: Phenotype data (sample metadata).
   - `CellVis_fullEnrichment_2023_12_16.csv`: GO-based subcellular annotations for ths dataset.

2. **RDS Files** (if available for reproducibility):
   - `paramssvm.rds`
   - `svmres.rds`
   - `paramsrf.rds`
   - `rfres.rds`
   - `paramsxgb.rds`
   - `xgbres.rds`
   - `inter.rds`
   - `my_csv_string_svm.rds`
   - `my_csv_string_rf.rds`
   - `my_csv_string_xgb.rds`
  

4. **R Libraries**:
   - `dplyr`, `tidyr`, `stringr`, `ggplot2`, `MSnbase`, `pRoloc`, `STRINGdb`, `subcellularvis`, `readr`, `purrr`.

## Script Workflow
The script consists of the following key steps:

### 1. **Data Preparation**
- **Input Files**: Reads the expression data (`f1.csv`), feature metadata (`f2.csv`), and phenotype data (`f3.csv`) to create the `MSnSet` object (`Control_trial`).
- **Validation**: Ensures column names in the expression data match sample names in the phenotype data.

### 2. **Machine Learning Predictions**
- Runs three ML algorithms:
  - **Support Vector Machine (SVM)**
  - **Random Forest (RF)**
  - **XGBoost (XGB)**
- Outputs:
  - Combined ML predictions for all proteins with their respective scores.

### 3. **STRINGdb Protein-Protein Interaction Analysis**
- Maps protein identifiers to STRINGdb IDs.
- Extracts PPI data and calculates combined interaction scores.
- Outputs:
  - PPI data combined with ML predictions and STRINGdb-derived FDR cutoffs.

### 4. **GO-Based Subcellular Annotations**
- Imports GO annotations from `CellVis_fullEnrichment_2023_12_16.csv`.
- Merges these annotations with protein predictions.
- Summarizes category-specific subcellular annotations.
- Outputs:
  - GO-based annotations for all proteins and GOCC-derived FDR cutoffs.

### 5. **Category-Specific FDR Cutoffs**
- Calculates **False Discovery Rate (FDR)** cutoffs for each category and ML algorithm:
  - STRINGdb-based FDR for SVM, RF, and XGB.
  - GO-based FDR for SVM, RF, and XGB.
- Combines STRINGdb and GO annotations to derive final category-specific score cutoffs for SVM, RF and XGB.


## Notes for Reproducibility
- Ensure all input files (`f1.csv`, `f2.csv`, `f3.csv`, `CellVis_fullEnrichment_2023_12_16.csv`) are correctly formatted and in the working directory.
- Column names in `f1.csv` must match sample names in `f3.csv`.
- RDS files can be regenerated if missing by rerunning the script.

