# CST Lung Cancer Viz

## Overview

This repo contains

* Lung cancer cell line PTM and gene expression data from CST and CCLE, respectively
* Data-Processing Python scripts
* the page [CST_Lung_Cancer_Viz](http://maayanlab.net/CST_Lung_Cancer_Viz/) source code

## Lung Cancer Data
Our collaborators at Cell Signaling Technology ([CST](https://www.cellsignal.com/)) used SILAC mass spectrometry to measure differential phosphorylation, acetylation, and methylation in a panel of 42 lung cancer cell lines compared to non-cancerous lung tissue. Gene expression data from 37 of these lung cancer cell lines was obtained from the publically available Cancer Cell Line Encyclopedia ([CCLE](https://portals.broadinstitute.org/ccle/home)).

## Data Processing Scripts
All data was processed using Python scripts in two broad steps: 1) data was pre-processed (e.g. calculating ratios of cancer vs non-cancer levels) and combined into a simple tab-separated format and 2) data was normalized and filtered in order to make heatmap visualizations.

Visualizations were made using the [Clustergrammer](https://github.com/MaayanLab/clustergrammer) web-based visualization tool. The [clustergrammer](clustergrammer) python module was used to normalize/filter data and produce JSONs for [clustergramer.js](js/clustergrammer.js) .

## Webpage Source Code
This repo contains the souce code for the site [CST_Lung_Cancer_Viz](http://maayanlab.net/CST_Lung_Cancer_Viz/). The page can also be seen through github.io [https://maayanlab.github.io/CST_Lung_Cancer_Viz/](https://maayanlab.github.io/CST_Lung_Cancer_Viz/) with some limited capabilities (due to github.io HTTPS requirements).

