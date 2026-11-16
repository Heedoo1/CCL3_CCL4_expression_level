# Figure 7c–g: data and plotting code

## Overview

This repository creates the data tables and individual plots for Figure 7c–g using GEO studies [GSE145926](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145926) and [GSE58682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58682). It produces separate panel files; assembling the full multipanel Figure 7 is not automated.

The plotting entry point uses the supplied manuscript plotting inputs. Figure 7e and 7f,g use donor/sample-level statistical tests. Figure 7d CCL3 and CCL4 use log-normalized expression in both Control and ALI. The Control colorbar displays 1–2 and the ALI colorbar 1–5; positive gene values outside each colorbar range are saturated for display only and are not changed or excluded; zero-expression cells are shown in grey. The same underlying log-normalization metric is used, but the same color does not represent the same expression value across the two panels because the displayed limits differ. The co-expression panels use cohort-specific 0–1 product scores and are relative, not directly comparable in absolute magnitude between cohorts.

## Requirements

- R 4.5.2 and the packages listed in `requirements-R.txt`. The listed versions were used on Windows; exact checks in the GSE145926 data-generation code can reject other versions. `ragg` is optional; if installed it is used for ALI PNG rendering.
- Python 3.10 or later with NumPy and SciPy from `requirements.txt` for GSE58682 data generation.
- Internet access and sufficient disk space for installing dependencies and downloading GEO inputs. The GSE145926 analyses require substantial RAM and time. Other operating systems and numerical libraries can yield different UMAP coordinates even with the same seed.

## Input data

The following files are included and are the default inputs to `00_RUN_PLOTS.R`:

| File | Purpose |
|---|---|
| `data/manuscript_plot_inputs/Figure7cd_Control_Manuscript.csv` | 20,327 Control cells for Figure 7c,d |
| `data/manuscript_plot_inputs/Figure7cd_ALI_Manuscript.csv` | 23,918 ALI cells and the submitted ALI UMAP coordinates for Figure 7c,d |
| `data/reference_workbooks/Figure7e_Source_Data.xlsx` | 12,937 independently selected AM cells, donor medians and Welch/BH statistics |
| `data/reference_workbooks/Figure_7f_7g_plotting_data_GSE58682.xlsx` | 56 exposure profiles and Figure 7f,g values |
| `data/reference_workbooks/Supplementary_Data_6.xlsx` | The associated Figure 7c,d per-cell supplementary data |

To regenerate data, download 12 original GEO H5 count matrices named `GSM4475048_C51_filtered_feature_bc_matrix.h5`, `GSM4475049_C52_filtered_feature_bc_matrix.h5`, `GSM4475050_C100_filtered_feature_bc_matrix.h5`, `GSM4339769_C141_filtered_feature_bc_matrix.h5`, `GSM4339770_C142_filtered_feature_bc_matrix.h5`, `GSM4339771_C143_filtered_feature_bc_matrix.h5`, `GSM4339772_C144_filtered_feature_bc_matrix.h5`, `GSM4339773_C145_filtered_feature_bc_matrix.h5`, `GSM4339774_C146_filtered_feature_bc_matrix.h5`, `GSM4475051_C148_filtered_feature_bc_matrix.h5`, `GSM4475052_C149_filtered_feature_bc_matrix.h5`, and `GSM4475053_C152_filtered_feature_bc_matrix.h5`. Place them together in an input directory. Download `GSE58682_series_matrix.txt.gz` for GSE58682 separately; the latter is the processed GEO Series Matrix, not raw CEL files. Neither download is included.

The GEO-generated Figure 7c,d workbook has two sheets, `Control` and `ALI`, with barcode, sample, cell type, UMAP coordinates, CCL3/CCL4 expression and 0–1 scores. The default plots instead read the included manuscript CSVs. The GEO workflow computes its own ALI v1.3 UMAP; its coordinates differ from the supplied manuscript ALI CSV. Data generation does not overwrite or replace the manuscript plotting inputs. The Control manuscript plot applies one common scale and translation to both axes for display. For the included input: `display_x = 0.7402556035793831 * UMAP_1 + 1.7165796921799406` and `display_y = 0.7402556035793831 * UMAP_2 + 0.7933126412823628`. These are display coordinates, so the manuscript plot axis numbers are not the original Source Data coordinates. The script derives the transform from the input ranges and fixed plot frame; the CSV itself is unchanged. An optional separate output plots the original coordinates exactly.

The Control source classifications are Alveolar macrophage-like (19,263), T cell (1,038), B cell-like (14), and Unresolved (12); there are **zero** Control neutrophils. To preserve the submitted Figure's title and legend, the default plot uses historical display labels: the label `Alveolar macrophage` denotes the source's macrophage-like cells, `B cell` denotes B cell-like cells, `Neutrophil` is a legend entry with no Control cells, and the 12 Unresolved cells are drawn in neutral grey without an entry in that historical legend. Display labels do not change, create, or reassign cells. The optional original-coordinate plot instead displays the source classifications in its legend.

## Installation and execution

1. Extract the ZIP and open an R console in the extracted repository directory. Install the pinned R packages and dependencies with:

   ```r
   install.packages("pak")
   packages <- scan("requirements-R.txt", what = "character", quiet = TRUE)
   pak::pkg_install(packages)
   ```

   Installing archived package versions may require compilation tools such as Rtools on Windows. Existing installations can be used if the versions match.

2. In a terminal opened in the repository directory, install the Python requirements for GSE58682 data generation:

   ```sh
   python -m pip install -r requirements.txt
   ```

3. To generate the submitted-input plots, run:

   ```sh
   Rscript 00_RUN_PLOTS.R
   ```

   To choose another empty output directory, append its path. Use a new empty output folder for each run; existing output is not overwritten by the launcher. To also make separate Control plots with the original input coordinates and source classification legend:

   ```sh
   Rscript 00_RUN_PLOTS.R /path/to/new_output --include-control-source
   ```

   To generate only the original-coordinate Control plots, use `Rscript 00_RUN_PLOTS.R /path/to/new_output --control-source-only`. These graphs are saved under `Control_original_coordinates`; the manuscript-style graphs are saved separately under `Control`.

4. Optional: regenerate the GSE145926 data from the 12 H5 files, using a new empty directory outside the input directory:

   ```sh
   Rscript code/data_generation/00_RUN_ALL.R /path/to/GSE145926_H5 /path/to/new_output
   ```

   The output `Figure7e_Source_Data.xlsx` has the same sheet contract as the included Figure 7e plotting input. To plot that file independently, set `FIG7E_INPUT_XLSX` to its path and `FIG7_OUTPUT_DIR` to a new output directory before running `Rscript code/plotting/06_Plot_Fig7e.R`. The default launcher continues to use the included manuscript inputs.

5. Optional: regenerate GSE58682 data from its matrix into a new folder:

   ```sh
   python code/data_generation/04_GSE58682_to_Excel.py /path/to/GSE58682_series_matrix.txt.gz /path/to/new_output
   ```

   Its `Figure_7f_7g_plotting_data_GSE58682.xlsx` has the same sheet contract as the included 7f,g input. To plot it independently, set `FIG7_INPUT_XLSX` and `FIG7_OUTPUT_DIR` before running `Rscript code/plotting/07_Plot_Fig7fg.R`.

Quote paths containing spaces. These commands read GEO files in place; they do not download or modify the inputs. No previously supplied Excel is used to calculate new GEO data.

## Output files

- Plotting writes PDF and 600-dpi PNG files below `figures/Control`, `figures/ALI`, `figures/Fig7e` and `figures/Fig7fg`, plus plot statistics and point tables. With a custom output directory these folders are created there.
- GSE145926 generation writes `Supplementary_Data_6.xlsx` (`per_cell_H.C`, `per_cell_ALI`, six columns each), `Figure7cd_Source_Data.xlsx` (`Control`, `ALI`), and `Figure7e_Source_Data.xlsx` (`Fig7e_Cell`, `Fig7e_Sample_Medians`, `Fig7e_Welch_BH`). Supplementary Data 6 contains the displayed 7c,d cells; Figure 7e uses a separate AM selection.
- GSE58682 generation writes `Supplementary_Data_7_GSE58682_reproduced.xlsx`, `Source_Data_Figure_7fg_extract.xlsx`, and `Figure_7f_7g_plotting_data_GSE58682.xlsx`.
- Figure 7e significance marks come from two-sided Welch tests of 3 Control versus 9 ALI donor medians, with BH correction across three tests. Figure 7f shows 18 Clean Air and 38 Ozone profiles: CCL3 and CCL4 marks use individual two-sided Welch tests, and the third mark uses a two-response CCL3/CCL4 Pillai MANOVA, although its box displays the mean of the two gene z-scores. Figure 7g reports two-sided Pearson correlations. The activation score includes CCL3 and CCL4 among its component genes.
