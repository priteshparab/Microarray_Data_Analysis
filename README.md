# Microarray_Data_Analysis
Microarray data analysis using R

## Introduction

Microarray technology is a high-throughput method used to measure the expression levels of thousands of genes simultaneously. It allows researchers to compare gene expression patterns between different conditions, such as treated vs. control samples, disease vs. healthy tissues, or different time points.

This repository provides an R-based workflow for analyzing microarray data, including:

- **Data normalization** using Robust Multi-array Average (RMA)
- **Visualization** of raw and normalized data
- **Statistical testing** (t-tests) to identify significant gene expression changes
- **Calculation of fold changes** for differential gene expression analysis

## Dataset Overview

**GSE21138**: Gene expression profiles in BA46 of subjects with schizophrenia and matched controls.

- **Organism:** Homo sapiens  
- **Experiment type:** Expression profiling by array  
- **Samples:** 30 subjects with schizophrenia. 

## Data Download and Preparation

- **Data Download:**  
  - Downloaded a total of **6 .CEL files** from the GEO database: 2 controls and 4 treateD samples.  
  - Also downloaded the **GPL annotation file** to map probe IDs to gene expression values.  

- **Data Organization:**  
  - Extracted the `.CEL` files and renamed them systematically for clarity and ease of use.  
  - Placed all `.CEL` files in a dedicated folder for analysis.  

- **R Setup and Analysis:**  
  - Set the working directory in R to the folder containing all `.CEL` files.  
  - Ran the **`Microarray_Script.R`** to perform data normalization, statistical analysis, and export results.

## Output and Filtering

- **Output Files:**  
  - Running the `Microarray_Script.R` generates an Excel file named **`microarray.xlsx`**, which contains:  
    - Normalized gene expression values  
    - p-values from t-tests  
    - Calculated fold changes for treated vs. control samples  

- **Filtering Significant Genes:**  
  - To identify differentially expressed genes, apply the following cutoffs:  
    - **Log2 Fold Change (|Log2FC|)**: ≥ 1 (upregulated or downregulated)  
    - **p-value**: ≤ 0.05  



