# CST Lung Cancer Viz

This repo contains
* Lung cancer cell line PTM and gene expression data from CST and CCLE, respectively
* Data-Processing Python scripts
* IPython/Jupyter [notebooks](notebooks) describing processing/normalization
* the source code for the page [CST_Lung_Cancer_Viz](http://maayanlab.net/CST_Lung_Cancer_Viz/)

# Lung Cancer Data
Our collaborators at Cell Signaling Technology ([CST](https://www.cellsignal.com/)) used SILAC mass spectrometry to measure differential phosphorylation, acetylation, and methylation in a panel of 42 lung cancer cell lines compared to non-cancerous lung tissue (the primary data is [here](lung_cellline_3_1_16)). Gene expression data from 37 of these lung cancer cell lines was obtained from the publically available Cancer Cell Line Encyclopedia ([(CCLE)](https://portals.broadinstitute.org/ccle/home)).

# Overview of CST PTM Data
The PTM data from CST is overviewed in two stages: 1) overview of the data, 2) overview of the normalization procedure

## Overview of the Data and Missing Data
The notebook, [CST_PTM_Data_Overview.ipynb](notebooks/CST_PTM_Data_Overview.ipynb), overviews the PTM data and discusses the reasoning behind our method of normanlization and filtering.

## Overview of Normalization
The notebook, [CST_PTM_Normalization_Overiew.ipynb](notebooks/CST_PTM_Normalization_Overview.ipynb)

The notebook, [t-SNE_Cell_Line_Clustering_Phosphorylation.ipynb](t-SNE_Cell_Line_Clustering_Phosphorylation.ipynb), uses the dimensionality reduction algorithm, t-SNE, to visualize cell line clustering based on phosphorylation data before and after normalization and filtering.

# Overview of CCLE Gene Expression Data
The python notebook, [CCLE_Gene_Expression_Data_Overview](notebooks/CCLE_Gene_Expression_Data_Overview.ipynb), overviews the gene expression data and discusses some results from our analysis of this data.

The notebook, [t-SNE_Cell_Line_Clustering_Gene_Expression.ipynb](notebooks/t-SNE_Cell_Line_Clustering_Gene_Expression.ipynb), uses t-SNE to visualize cell line clustering based on gene expression data.

# Data Processing Scripts
All data was processed using Python scripts in two broad steps: 1) data was pre-processed (e.g. calculating ratios of cancer vs non-cancer levels) and combined into a simple tab-separated format and 2) data was normalized and filtered in order to make heatmap visualizations. Visualizations for the webpage were made using the [Clustergrammer](https://github.com/MaayanLab/clustergrammer) web-based visualization tool. The [clustergrammer](clustergrammer) python module was used to normalize/filter data and produce JSONs for [clustergramer.js](js/clustergrammer.js) .

## Data Pre-processing
The script [process_latest_cst_data.py](process_latest_cst_data.py) was used to 1) calculate the log2 peak ratios of lung cancer cell lines to their associated Normal Pool from the same plex, and 2) combine these ratios into a single tab-separated file.

## Data Normalization and Filtering
The script [make_cst_homepage_figures.py](make_cst_homepage_figures.py) was used to make Clustergrammer interactive visualizations for the website. Specifically, this script normalizes and processes the data and creates the [JSONs](json) for the front-end visualizations.


# Webpage Source Code
This repo contains the souce code for the site [CST_Lung_Cancer_Viz](http://maayanlab.net/CST_Lung_Cancer_Viz/). The page can also be seen through github.io [https://maayanlab.github.io/CST_Lung_Cancer_Viz/](https://maayanlab.github.io/CST_Lung_Cancer_Viz/) with some limited capabilities (due to github.io HTTPS requirements).

