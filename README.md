# Code and data for the manuscript

This repository provides R code and data to reproduce results and figures from the manuscript:

> Keck, F. et al. Extracting massive ecological data on state and interactions of species using large language models (2025).

## System requirements

The code was developed and tested on Linux. Some tools used are specific to this platform. The details about the R environment on which the code has been tested is fully described in `/session_info.txt`.

## Installation

Install R packages from CRAN and GitHub:

```         
install.packages(c("tidyverse", "openalexR", "httr2", "dplyr",
"readr", "tidyr", "curl", "jsonlite", "stringdist", "cli",
"R.utils", "tidygraph", "ggraph", "patchwork", "tidyheatmaps",
"RColorBrewer", "pheatmap", "igraph", "htmltools",
"xml2", "magrittr", "taxize", "pluralize", "stringr", "wikitaxa",
"rvest"))

remotes::install_github("fkeck/flexitarian")
```

Using the package's binaries provided by CRAN, installation takes a few minutes on standard computer. The versions of the packages used to generate the results of the manuscript are provided in `/session_info.txt`.

`Python 3.10+` and the libraries `spaCy` and `TaxoNERD` with the model `en_core_eco_biobert_weak` are needed to perform NER.

Linux commands `tar` and `jq` are also required.

## Instructions

### Data

Data can be found in the `/data` directory. The main data file (data/`save_R_gdata_2.csv.tar.gz`) must be uncompressed to run the analyses. The raw text content of the processed publications can be found on [PubMed OA](https://pmc.ncbi.nlm.nih.gov/tools/openftlist/). The scripts to download, preprocess and generate the raw and intermediate files are in the `R/` directory.

### Analyses

Scripts can be found in the `/R` directory. They are organized in separated numbered modules. The complete execution can take several hours. The code does not include the interactions with the OpenAI API platform which is subject to a charge.
