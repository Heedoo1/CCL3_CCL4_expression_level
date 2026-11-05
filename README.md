# CCL3_CCL4_expression_level

## Overview
This repository contains the bioinformatics analysis code used in the study  
“Aging Reconfigures Extracellular Vesicle Signaling Driving Environmental Inflammatory Susceptibility”.

All bioinformatics analyses represented in this repository were performed using the R programming environment.

This repository contains analysis scripts only and does not include wet-lab experimental procedures, sample preparation protocols, or primary data generation.
The R scripts document the computational analyses performed for the manuscript, including UMAP-based cell-type visualization and chemokine expression comparisons under inflammatory and environmental exposure conditions.

## Data Sources
The bioinformatics analyses in this repository are based on publicly available datasets. 
Transcriptomic data were obtained from the Gene Expression Omnibus (GEO). GEO accession numbers: GSE145926, GSE58682. 
These analyses represent secondary bioinformatics analyses of existing datasets. No new experimental data were generated as part of the computational analyses documented in this repository.

## System Requirements

### Operating system
- Windows 11 x64  (tested operating system; compatibility with other operating systems has not been independently verified)

### R version
- R ≥ 4.2.0 (tested on R 4.5.2)

### Required R packages
- Seurat
- dplyr
- ggplot2
- patchwork
- scales
- BiocManager
- oligo
- pd.hugene.1.0.st.v1
- hugene10sttranscriptcluster.db
- AnnotationDbi
- stringr
- writexl
- tidyr

## Installation Guide 

1. Install R (version ≥ 4.2.0) from CRAN.
2. Install the required R packages:

```r
install.packages(c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "scales",
  "BiocManager", "stringr", "writexl", "tidyr"
))

BiocManager::install(c(
  "oligo",
  "pd.hugene.1.0.st.v1",
  "hugene10sttranscriptcluster.db",
  "AnnotationDbi"
))
```
## Instructions for Use
The R scripts in the `scripts` directory are organized by analysis component and numbered for reference.

1. 01_Control_UMAP_visualization.R
2. 02_ALI_UMAP_visualization.R
3. 03_CCL3_CCL4_expression_comparison.R
4. 04_Ozone_inhalation_CCL3_CCL4_expression_comparison.R
5. 05_Ozone_inhalation_CCL3_CCL4_activationScore_correlation.R

Each script documents the analysis indicated by its filename. Execution requires the corresponding input data and, where applicable, intermediate analysis objects generated during the original workflow.

## Analysis Inputs and Reproducibility

Some R scripts require intermediate data objects generated during the original analysis, including processed Seurat objects and analysis-specific data frames.

For example, `obj` and `am_cells` are used as input objects in the corresponding scripts.

These intermediate objects are not included in the current repository and must be provided before executing the relevant scripts.
