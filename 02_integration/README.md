# Beta-cell and lineage-specific reintegration

These plain R scripts reproduce the original object-construction workflow. They retain the sequential analysis style of the existing GitHub figure scripts. Run them in numeric order **from the same working directory** so relative input/output paths resolve.

| Script | Input | Principal saved output | Leiden resolution | Slingshot root |
|---|---|---|---|---|
| `01_beta_atlas_subset_integration_trajectory.R` | `scRNA_integrated_harmony_classified.rds` | `beta.integrated.rds`, `cds_all.rds`, `Beta_T2D_lineage1.rds`, `Beta_T2D_lineage2.rds`, `Beta_T2D_lineage3.rds` | 0.00004 | 3 |
| `02_lin1_reintegration.R` | `Beta_T2D_lineage1.rds` | `Lin1.beta.reintegrated.rds`, `cds_lin1.rds`, `slingshot_sim_lin1.rds`, `sce_lin1.rds` | 0.00004 | 2 |
| `03_lin2_reintegration.R` | `Beta_T2D_lineage2.rds` | `lin2.beta.reintegrated.rds`, `cds_lin2.rds`, `slingshot_sim_lin2.rds`, `sce_lin2.rds` | 0.00011 | 1 |
| `04_lin3_reintegration.R` | `Beta_T2D_lineage3.rds` | `lin3.beta.reintegrated.rds`, `cds_lin3.rds`, `slingshot_sim_lin3.rds` | 0.0001 | 1 |
| `05_alpha_atlas_subset_integration.R` | `scRNA_integrated_harmony_classified.rds` | `alpha.integrated.rds`, `cds_alpha.rds` | 0.00001 | none (no trajectory) |
| `06_tradeSeq_association_tests.R` | `slingshot_sim_lin1/2/3.rds`, `Lineage1/2/3_sling_perc_expr.tsv` | `sce_lin1_1/1_2/2_1/2_2/3.rds`, `*_pseudotime_start_end_association_filt.tsv` | — | — |
| `07_pseudobulk_counts.R` | `Lin1/lin2/lin3.beta.reintegrated.rds` | `*_pseudobulk_cts.tsv`, `*_pseudobulk_coldata.tsv` in `Figure3_S5/`, `Figure4_S6/`, `Figure5_S8/`, `Figure6_S10/` | — | — |

Original SCTransform, dimensional reduction, Monocle3, Slingshot and lineage subsetting logic is retained. The original R Markdown chunks marked `eval=FALSE` are now executable code, including object generation and `saveRDS()` calls. Optional exploratory plots are retained as in the originals. Original datasets and package dependencies are required; the outputs have not been numerically regenerated or compared to the archived manuscript objects.

Small necessary fixes are annotated with `#!`: explicitly choose `assigned_label` before the initial beta subset, select the `RNA` assay for the SCT normalization in `01` only (the lineage inputs in `02`-`04` keep their saved default assay, SCT, as in the original run), remove a trailing `RunUMAP` argument comma and correct two erroneous `mean_value2`/`median_value2` column references in lineage 3. No analysis thresholds or root assignments were changed. `cds_all.rds` is saved by `01` (the original script did not save it); `sce_all.rds` (tradeSeq fitGAM output read by Figure 2) is generated elsewhere.

`05` keeps datasets with more than 200 alpha cells (the beta-cell scripts use more than 100), as in the original run. `06` consolidates the five original per-branch tradeSeq scripts (fitGAM with study as covariate, nknots = 9; associationTest; startVsEndTest filtered to FDR <= 0.05 and expression in >= 5% of lineage cells). The joint two-branch fits `sce_lin1.rds` and `sce_lin2.rds` (read by the Figure 3, 4 and 5 scripts) are generated in `02` and `03`; those lines were reconstructed and verified against the archived objects (dimensions, both lineages, 9 knots, study covariate).

`07` sums raw RNA counts per study x condition x cluster (`monocle_clusters2`); each pseudobulk sample therefore pools all donors of one study and condition, and the number of cells per sample is stored in `counts`.
