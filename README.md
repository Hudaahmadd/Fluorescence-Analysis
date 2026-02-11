# Fluorescence Analysis Pipeline

## Overview
This repository contains an R-based analysis pipeline for processing fluorescence data from plate reader experiments.

The pipeline:
- Maps fluorescence values to strain/OD levels using a plate map (Plate_info.txt)
- Generates fluorescence distribution plots
- Performs Tukey HSD post-hoc statistical testing (adjusted p-values)
- Exports high-resolution TIFF figures and CSV tables

## Repository Structure
- plate_reader_analysis.r
- LICENSE
- data/
  - Plate_info/Plate_info.txt
  - plate_reader/fluorescence.csv
- outputs/ (generated automatically)

## System Requirements

### Tested Environment
- Operating System: macOS 15.7.3 (Build 24G419)
- R Version: 4.5.0 (2025-04-11)

The script is expected to run on Windows and Linux systems with compatible R installations.

### Required R Packages
- tidyverse
- ggplot2
- ggbeeswarm
- rstatix
- ggpubr

## Installation

### 1. Install R
Download R (>= 4.5.0 recommended) from:
https://cran.r-project.org/

### 2. Install Required Packages
Run in R:

    install.packages(c(
      "tidyverse",
      "ggplot2",
      "ggbeeswarm",
      "rstatix",
      "ggpubr"
    ))

Typical installation time: < 5 minutes on a standard desktop computer.

## Demo

### Running the Included Dataset
From the repository root, run:

    Rscript plate_reader_analysis.r

### Expected Output
The script generates:

- outputs/plate_reader_data/RUN_fluorescence/fluorescence_plots/
  - Flu_<Condition>.tiff
- outputs/plate_reader_data/RUN_fluorescence/fluorescence_plots_with_stats/
  - Flu_withStats_<Condition>.tiff
- outputs/plate_reader_data/RUN_fluorescence/CSV_outputs/
  - fluorescence_data.csv
  - TukeyHSD_Flu_<Condition>.csv

### Expected Runtime
~1.17 seconds total on a standard desktop computer (measured on macOS 15.7.3 with R 4.5.0 using the dataset included in this repository).

## Instructions for Use on New Data

1) Replace the input files:
- data/plate_reader/fluorescence.csv
- data/Plate_info/Plate_info.txt

2) Required columns

Plate_info.txt must contain:
- row
- column
- strain

fluorescence.csv must contain:
- Row (or row)
- Column (or column)
- fluorescence (or intensity)
- Optional: Condition

3) Run:

    Rscript plate_reader_analysis.r

## Statistical Analysis
- Fluorescence values are grouped by strain/OD.
- Tukey’s Honest Significant Difference (HSD) post-hoc testing is applied.
- Only adjusted p-values (p.adj < 0.05) are annotated on statistical plots.
- Full Tukey test results are exported as CSV files for transparency.

## License
This software is distributed under the MIT License (see LICENSE).
