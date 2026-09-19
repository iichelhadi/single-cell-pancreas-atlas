# Human pancreatic beta-cell atlas analysis

This repository contains the analysis scripts used to generate the computational results and figures for the revised manuscript describing an integrated human pancreatic single-cell atlas and beta-cell states across aging and diabetes.

## Figure mapping

- **Figure 1, Figures S1–S2** → `Figure_1_S1_S2/`
- **Figure 2, Figure S3** → `Figure2_S3/`
- **Figure 3, Figure S5** → `Figure3_S5/`
- **Figure 4, Figure S6** → `Figure4_S6/`
- **Figure 5, Figure S8** → `Figure5_S8/`
- **Figure 6, Figure S10** → `Figure6_S10/`
- **Figure 7, Figure S11** → `Figure7_S11/`

Some revised analyses are implemented in dedicated scripts within the corresponding figure directories, including tradeSeq trajectory-module analyses, donor-level statistical analyses, sensitivity analyses and SCENIC consensus analyses.

## Data

Processed single-cell objects and computational source data underlying the manuscript figures are available from Zenodo:

DOI: 10.5281/zenodo.22767582

The public scRNA-seq datasets used in the study are available from GEO and ENA under the accession numbers reported in the manuscript.

## Code

Scripts use relative paths where possible and are organized by manuscript figure.

The mutual information network analysis used in Figure 7 is implemented using the `miEdgeR` package:

https://github.com/iichelhadi/miEdgeR

## Reproducibility

Random seeds are specified for stochastic analyses where applicable. Package and software versions used for the analyses are reported in the manuscript Methods and associated computational environment information.

© 2026 Elhadi Iich
