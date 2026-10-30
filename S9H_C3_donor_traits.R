################################################################################
### S9H: Lineage 3 cluster C3 donor-level trait associations
###
### CONSISTENCY AUDIT ITEM (not an explicit reviewer request).
### R1-9 objected to the donor-trait statistical methodology for C2/C5. The
### identical analysis exists for Lineage 3 C3 in Figure6_FigureS8_clean.R
### (section 9, "Figure S8H" in the script / S9H in the manuscript) and carries
### the same flaws, in a more severe form:
###
###   Original S9H                          This script
###   ------------------------------------  ------------------------------------
###   donor unit = raw orig.ident           harmonized donor_id (annotation
###     (each pooled GSE = ONE "donor")       joins; GSE81608 donor_id,
###                                           GSE153855 BioSample, metadata
###                                           fallback for GSE98887/GSE124742)
###   Pearson cor.test on proportions       Spearman (proportions are bounded,
###                                           non-normal)
###   t.test for sex                        Wilcoxon rank-sum
###   no multiple-testing correction        Benjamini-Hochberg
###   stats_table hard-coded                computed from the live result
###
### Mirrors R1-9__R2-15.R exactly except for the lineage, cluster, and
### BH family. Do not diverge from that script's methodology here — the whole
### point is that the two tables are held to one standard.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
})

set.seed(123)

dir.create("revision_out", showWarnings = FALSE)

################################################################################
### 0. Configuration
################################################################################

#! BH FAMILY DECISION — the one open methodological choice.
#!   "separate" (default) : correct across this script's 5 C3 tests only.
#!                          C2/C5 keeps its own 10-test correction. Each figure
#!                          answers a distinct question about a distinct
#!                          lineage; state the family explicitly in each legend.
#!   "pooled"             : you must instead re-run R1-9__R2-15.R and this
#!                          script together over all 15 tests. NOT supported by
#!                          running this script alone — it would silently
#!                          correct 5 tests while claiming 15.
BH_FAMILY <- "separate"
stopifnot(BH_FAMILY %in% c("separate", "pooled"))
if (BH_FAMILY == "pooled") {
  stop(
    "Pooled correction requires a combined run across C2/C5 and C3 together. ",
    "This script corrects within its own 5 tests only. Either set BH_FAMILY ",
    "to 'separate', or merge the two result tables and p.adjust() once ",
    "downstream — and then re-issue the C2/C5 numbers, which will change."
  )
}

LIN3_OBJ  <- "/data/scratch/gmsv2524/revision/Figure6_S8/lin3.beta.reintegrated.rds"  #! confirm path
CLUSTER_COL <- "monocle_clusters2"   #! confirmed from Figure6_FigureS8_clean.R
TARGET_CLUSTER <- "3"                #! C3, the pre-T2D cluster in Lineage 3

# Annotation CSVs must be present in ./revision_out/ (same files used by
# R1-9__R2-15.R). Copy them over if running from a different directory.

################################################################################
### 1. Load Lineage 3 object
################################################################################

lin3_obj <- readRDS(LIN3_OBJ)
print(dim(lin3_obj))

################################################################################
### 2. Extract current cell metadata
################################################################################

required_cols <- c("orig.ident", "sample_name", "Study", CLUSTER_COL,
                   "Age", "BMI", "Hba1c", "Gender", "Ethnicity", "Condition")
missing_cols <- setdiff(required_cols, colnames(lin3_obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Lineage 3 object is missing expected metadata columns: ",
       paste(missing_cols, collapse = ", "),
       ". Do not proceed — the Lineage 2 script assumed all of these.")
}

metadata <- lin3_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::transmute(
    cell_id,
    orig.ident = as.character(orig.ident),
    sample_name = as.character(sample_name),
    Study = as.character(Study),
    cluster = as.character(.data[[CLUSTER_COL]]),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI),
    Hba1c = as.numeric(Hba1c),
    Gender = as.character(Gender),
    Ethnicity = as.character(Ethnicity),
    Condition = as.character(Condition)
  )

stopifnot(
  nrow(metadata) == ncol(lin3_obj),
  !anyDuplicated(metadata$cell_id)
)

cat("\nCells per Lineage 3 cluster:\n")
print(table(metadata$cluster, useNA = "ifany"))

## --- Raw categorical value check ------------------------------------------
## Across 18 integrated studies, capitalization/whitespace in Gender and
## Ethnicity cannot be assumed consistent. This matters in TWO places:
##   1. the sex/ethnicity test filters (handled by normalizing inside the
##      test functions below), and
##   2. donor_signature, which defines donor groups for the fallback studies
##      GSE98887/GSE124742.
## Normalizing before building donor_signature would group fallback donors by
## a different rule than R1-9__R2-15.R used for C2/C5 — reintroducing the very
## inconsistency this rerun exists to remove. So signatures are built on RAW
## values (matching C2/C5), and normalization is confined to the tests.
## If the tables below show variant spellings, STOP: the C2/C5 run needs the
## same inspection before either result is final.
cat("\nRaw Gender values (check for variant spellings):\n")
print(table(metadata$Gender, useNA = "ifany"))
cat("\nRaw Ethnicity values (check for whitespace/spelling variants):\n")
print(table(metadata$Ethnicity, useNA = "ifany"))

if (!TARGET_CLUSTER %in% metadata$cluster) {
  stop("Cluster '", TARGET_CLUSTER, "' not found in ", CLUSTER_COL,
       ". Check which label denotes C3 in this object.")
}

################################################################################
### 3. Read source annotation mappings
################################################################################

ann81608 <- readr::read_csv("./revision_out/GSE81608.csv", show_col_types = FALSE) %>%
  dplyr::transmute(
    sample_name = as.character(Run),
    annotated_donor = as.character(donor_id)
  ) %>%
  dplyr::filter(!is.na(sample_name), sample_name != "",
                !is.na(annotated_donor), annotated_donor != "") %>%
  dplyr::distinct()

ann153855 <- readr::read_csv("./revision_out/GSE153855.csv", show_col_types = FALSE) %>%
  dplyr::transmute(
    sample_name = as.character(Run),
    annotated_biosample = as.character(BioSample)
  ) %>%
  dplyr::filter(!is.na(sample_name), sample_name != "",
                !is.na(annotated_biosample), annotated_biosample != "") %>%
  dplyr::distinct()

################################################################################
### 4. Verify that each Run has only one source identifier
################################################################################

stopifnot(
  nrow(ann81608 %>% dplyr::count(sample_name) %>% dplyr::filter(n > 1)) == 0,
  nrow(ann153855 %>% dplyr::count(sample_name) %>% dplyr::filter(n > 1)) == 0
)

################################################################################
### 5. Join source identifiers to current cells
################################################################################

metadata <- metadata %>%
  dplyr::left_join(ann81608, by = "sample_name") %>%
  dplyr::left_join(ann153855, by = "sample_name")

################################################################################
### 6. Verify annotation joins
###
### NOTE: unlike R1-9__R2-15.R, the expected donor counts are NOT hard-coded.
### Lineage 3 contains a different cell subset than Lineage 2, so the number of
### donors represented WILL differ from 16/11/2/4. Hard-coding the Lineage 2
### numbers here would either crash the script or pass for the wrong reason.
### Inspect the printed table and confirm the match rate is 100%.
################################################################################

join_check <- metadata %>%
  dplyr::filter(orig.ident %in% c("GSE81608", "GSE153855")) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_matched = dplyr::case_when(
      dplyr::first(orig.ident) == "GSE81608"  ~ sum(!is.na(annotated_donor)),
      dplyr::first(orig.ident) == "GSE153855" ~ sum(!is.na(annotated_biosample))
    ),
    n_identifiers = dplyr::case_when(
      dplyr::first(orig.ident) == "GSE81608"  ~ dplyr::n_distinct(annotated_donor, na.rm = TRUE),
      dplyr::first(orig.ident) == "GSE153855" ~ dplyr::n_distinct(annotated_biosample, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("\nAnnotation join check (verify n_matched == n_cells for each):\n")
print(join_check)

if (nrow(join_check) > 0 && any(join_check$n_cells != join_check$n_matched)) {
  stop("Incomplete annotation join — some cells did not receive a donor/BioSample ID.")
}

################################################################################
### 7. Construct fallback metadata signatures
################################################################################

fallback_studies <- c("GSE98887", "GSE124742")

signature_value <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x
}

metadata <- metadata %>%
  dplyr::mutate(
    donor_signature = paste(
      signature_value(Age), signature_value(BMI), signature_value(Hba1c),
      signature_value(Gender), signature_value(Ethnicity), signature_value(Condition),
      sep = "|"
    )
  )

################################################################################
### 8. Create final harmonized donor ID
################################################################################

metadata <- metadata %>%
  dplyr::mutate(
    donor_id = dplyr::case_when(
      orig.ident == "GSE81608"  ~ paste0("GSE81608::", annotated_donor),
      orig.ident == "GSE153855" ~ paste0("GSE153855::", annotated_biosample),
      orig.ident %in% fallback_studies ~ paste0(orig.ident, "::metadata::", donor_signature),
      TRUE ~ paste0(Study, "::", orig.ident)
    ),
    donor_id_source = dplyr::case_when(
      orig.ident == "GSE81608"  ~ "annotated donor ID",
      orig.ident == "GSE153855" ~ "annotated BioSample",
      orig.ident %in% fallback_studies ~ "metadata-defined group",
      TRUE ~ "existing donor/sample identifier"
    )
  )

stopifnot(!any(is.na(metadata$donor_id)), !any(metadata$donor_id == ""))

################################################################################
### 9. Inspect final donor counts (report, do not assert)
################################################################################

donor_counts_by_study <- metadata %>%
  dplyr::group_by(Study, orig.ident, donor_id_source) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_donor_groups = dplyr::n_distinct(donor_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Study, orig.ident)

cat("\nDonor counts for the four pooled studies (Lineage 3 subset):\n")
print(
  donor_counts_by_study %>%
    dplyr::filter(orig.ident %in% c("GSE81608", "GSE153855", "GSE98887", "GSE124742"))
)
cat("\n#! These will NOT match the Lineage 2 counts (16/11/2/4). Different cells.\n")

################################################################################
### 10. Verify metadata consistency within final donor IDs
################################################################################

metadata_inconsistency <- metadata %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(
    n_studies = dplyr::n_distinct(Study, na.rm = TRUE),
    n_age = dplyr::n_distinct(Age, na.rm = TRUE),
    n_BMI = dplyr::n_distinct(BMI, na.rm = TRUE),
    n_HbA1c = dplyr::n_distinct(Hba1c, na.rm = TRUE),
    n_Gender = dplyr::n_distinct(Gender, na.rm = TRUE),
    n_Ethnicity = dplyr::n_distinct(Ethnicity, na.rm = TRUE),
    n_Condition = dplyr::n_distinct(Condition, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(
    n_studies > 1 | n_age > 1 | n_BMI > 1 | n_HbA1c > 1 |
      n_Gender > 1 | n_Ethnicity > 1 | n_Condition > 1
  )

cat("\nDonor IDs with inconsistent metadata (should be empty):\n")
print(as.data.frame(metadata_inconsistency))
stopifnot(nrow(metadata_inconsistency) == 0)

################################################################################
### 11. Create one metadata row per donor/group
################################################################################

first_numeric <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[[1]]
}

first_character <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  as.character(x[[1]])
}

donor_meta <- metadata %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(
    Study = dplyr::first(Study),
    orig.ident = dplyr::first(orig.ident),
    donor_id_source = dplyr::first(donor_id_source),
    Age = first_numeric(Age),
    BMI = first_numeric(BMI),
    Hba1c = first_numeric(Hba1c),
    Gender = first_character(Gender),
    Ethnicity = first_character(Ethnicity),
    Condition = first_character(Condition),
    .groups = "drop"
  )

cat("\nTotal donor/groups in Lineage 3:\n")
print(dplyr::n_distinct(donor_meta$donor_id))

cat("\nDonor ID sources:\n")
print(table(donor_meta$donor_id_source))

# Trait coverage — mirrors the GSE153855 all-NA finding from the C2/C5 run.
# Studies with no usable covariates cannot contribute to any trait test, and
# that fact belongs in the rebuttal rather than being left implicit.
cat("\nTrait coverage by study (n donors with non-missing values):\n")
print(
  donor_meta %>%
    dplyr::group_by(orig.ident) %>%
    dplyr::summarise(
      n_donors = dplyr::n(),
      n_Age = sum(!is.na(Age)),
      n_BMI = sum(!is.na(BMI)),
      n_HbA1c = sum(!is.na(Hba1c)),
      n_Gender = sum(!is.na(Gender) & Gender != ""),
      n_Ethnicity = sum(!is.na(Ethnicity) & Ethnicity != ""),
      .groups = "drop"
    ),
  n = Inf
)

################################################################################
### 12. Calculate donor-level C3 occupancy
################################################################################

donor_totals <- metadata %>%
  dplyr::count(donor_id, name = "total_lineage3_cells")

donor_props <- metadata %>%
  dplyr::filter(cluster == TARGET_CLUSTER) %>%
  dplyr::count(donor_id, cluster, name = "n_cells") %>%
  tidyr::complete(
    donor_id = unique(metadata$donor_id),
    cluster = TARGET_CLUSTER,
    fill = list(n_cells = 0)
  ) %>%
  dplyr::left_join(donor_totals, by = "donor_id") %>%
  dplyr::mutate(proportion = n_cells / total_lineage3_cells) %>%
  dplyr::left_join(donor_meta, by = "donor_id") %>%
  dplyr::arrange(donor_id)

stopifnot(
  all(!is.na(donor_props$total_lineage3_cells)),
  all(donor_props$total_lineage3_cells > 0),
  all(donor_props$proportion >= 0),
  all(donor_props$proportion <= 1)
)

cat("\nDonor occupancy dimensions:\n")
print(dim(donor_props))

################################################################################
### 13. Statistical test functions (identical to R1-9__R2-15.R)
###
### The stopifnot guards from the C2/C5 script are replaced with graceful
### returns: Lineage 3 has fewer donors, and a trait with only one level or
### too few donors should be reported as untestable rather than crashing the
### run partway through.
################################################################################

continuous_test <- function(dat, variable) {
  d <- dat %>% dplyr::filter(!is.na(.data[[variable]]), !is.na(proportion))
  if (nrow(d) < 3) {
    return(tibble::tibble(
      trait = variable, test_method = "Spearman correlation",
      n_donors = nrow(d), statistic = NA_real_, statistic_name = "rho",
      p_value = NA_real_, group_summary = "insufficient donors"
    ))
  }
  fit <- stats::cor.test(d[[variable]], d$proportion,
                         method = "spearman", exact = FALSE)
  tibble::tibble(
    trait = variable, test_method = "Spearman correlation",
    n_donors = nrow(d), statistic = unname(fit$estimate),
    statistic_name = "rho", p_value = fit$p.value,
    group_summary = NA_character_
  )
}

sex_test <- function(dat) {
  # Normalize HERE ONLY, not upstream — donor_signature must stay on raw
  # values to match the C2/C5 donor grouping. See the diagnostic in section 2.
  d <- dat %>%
    dplyr::mutate(Gender = tolower(trimws(Gender))) %>%
    dplyr::filter(Gender %in% c("male", "female"), !is.na(proportion))
  if (dplyr::n_distinct(d$Gender) < 2 || nrow(d) < 3) {
    return(tibble::tibble(
      trait = "Gender", test_method = "Wilcoxon rank-sum",
      n_donors = nrow(d), statistic = NA_real_, statistic_name = "W",
      p_value = NA_real_, group_summary = "insufficient donors or single sex"
    ))
  }
  fit <- stats::wilcox.test(proportion ~ Gender, data = d, exact = FALSE)
  summary_by_group <- d %>%
    dplyr::group_by(Gender) %>%
    dplyr::summarise(n = dplyr::n(),
                     median_proportion = median(proportion, na.rm = TRUE),
                     .groups = "drop")
  tibble::tibble(
    trait = "Gender", test_method = "Wilcoxon rank-sum",
    n_donors = dplyr::n_distinct(d$donor_id), statistic = unname(fit$statistic),
    statistic_name = "W", p_value = fit$p.value,
    group_summary = paste0(summary_by_group$Gender, ": n=", summary_by_group$n,
                           ", median=", round(summary_by_group$median_proportion, 4),
                           collapse = "; ")
  )
}

ethnicity_test <- function(dat) {
  # Whitespace normalization only, same rationale as sex_test. Deliberately
  # NOT lowercased — ethnicity labels are proper nouns and case-folding could
  # merge categories the source studies intended as distinct.
  d <- dat %>%
    dplyr::mutate(Ethnicity = trimws(Ethnicity)) %>%
    dplyr::filter(!is.na(Ethnicity), Ethnicity != "", !is.na(proportion))
  if (dplyr::n_distinct(d$Ethnicity) < 2 || nrow(d) < 3) {
    return(tibble::tibble(
      trait = "Ethnicity", test_method = "Kruskal-Wallis",
      n_donors = nrow(d), statistic = NA_real_, statistic_name = "chi-square",
      p_value = NA_real_, group_summary = "insufficient donors or single group"
    ))
  }
  fit <- stats::kruskal.test(proportion ~ Ethnicity, data = d)
  tibble::tibble(
    trait = "Ethnicity", test_method = "Kruskal-Wallis",
    n_donors = dplyr::n_distinct(d$donor_id), statistic = unname(fit$statistic),
    statistic_name = "chi-square", p_value = fit$p.value,
    group_summary = NA_character_
  )
}

################################################################################
### 14. Run C3 association analyses
################################################################################

results <- dplyr::bind_rows(
  continuous_test(donor_props, "BMI"),
  continuous_test(donor_props, "Age"),
  continuous_test(donor_props, "Hba1c"),
  sex_test(donor_props),
  ethnicity_test(donor_props)
) %>%
  dplyr::mutate(cluster = "C3", .before = 1) %>%
  dplyr::mutate(p_adjusted_BH = stats::p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(match(trait, c("BMI", "Age", "Hba1c", "Gender", "Ethnicity")))

print(results, n = Inf)


################################################################################
### Revised Figure S9H table
################################################################################

stats_table_S9H <- results %>%
  dplyr::mutate(
    Metadata = dplyr::recode(
      trait,
      BMI = "BMI",
      Age = "Age",
      Hba1c = "HbA1c",
      Gender = "Sex",
      Ethnicity = "Ethnicity"
    ),
    Cluster = cluster,
    Test = test_method,
    Statistic = dplyr::case_when(
      statistic_name == "rho" ~ paste0(
        "rho=", sprintf("%.3f", statistic)
      ),
      statistic_name == "W" ~ paste0(
        "W=", sprintf("%.0f", statistic)
      ),
      statistic_name == "chi-square" ~ paste0(
        "chi²=", sprintf("%.3f", statistic)
      )
    ),
    P_value = sprintf("%.4f", p_value),
    P_adjusted_BH = sprintf("%.4f", p_adjusted_BH)
  ) %>%
  dplyr::select(
    Metadata,
    Cluster,
    Test,
    n_donors,
    Statistic,
    P_value,
    P_adjusted_BH
  )

stats_table_S9H

library(flextable)
library(officer)

figS9H_table <- stats_table_S9H %>%
  flextable::flextable() %>%
  flextable::set_header_labels(
    Metadata = "Trait",
    Cluster = "Cluster",
    Test = "Test",
    n_donors = "n donors",
    Statistic = "Statistic",
    P_value = "P",
    P_adjusted_BH = "BH-adjusted P"
  ) %>%
  flextable::autofit() %>%
  flextable::bg(
    i = seq(1, nrow(stats_table_S9H), by = 2),
    bg = "#FFFFFF",
    part = "body"
  ) %>%
  flextable::bg(
    i = seq(2, nrow(stats_table_S9H), by = 2),
    bg = "#f0f0f0",
    part = "body"
  ) %>%
  flextable::bg(
    bg = "#f0f0f0",
    part = "header"
  ) %>%
  flextable::bold(
    bold = TRUE,
    part = "header"
  ) %>%
  flextable::align(
    align = "center",
    part = "all"
  ) %>%
  flextable::hline(
    part = "all",
    border = officer::fp_border(
      color = "black",
      width = 0.5
    )
  )

figS9H_table


cat("\n--- Comparison with the superseded S9H values ---\n")
cat("Original (Pearson/t-test, orig.ident as donor, uncorrected):\n")
cat("  BMI       r=0.375, P=0.002\n")
cat("  Age       r=0.055, P=0.644\n")
cat("  HbA1c     r=0.191, P=0.170\n")
cat("  Gender    t=-0.500, P=0.619\n")
cat("  Ethnicity chi2=7.910, P=0.019\n")
cat("\nIf BMI and/or ethnicity no longer survive BH correction, the manuscript\n")
cat("claim that obesity and ethnic background contribute to C3 prevalence must\n")
cat("be removed, and the R1-10 lipotoxicity response revised accordingly.\n")
cat("If they DO survive: retain as an association only. Surviving significance\n")
cat("is not surviving causality — the 'reinforces lipotoxicity' gloss goes\n")
cat("either way.\n")

################################################################################
### 15. Save outputs
################################################################################

readr::write_tsv(results,      "revision_out/S9H_C3_donor_trait_associations.tsv")
readr::write_tsv(donor_props,  "revision_out/S9H_C3_donor_level_occupancy.tsv")
readr::write_tsv(donor_meta,   "revision_out/S9H_C3_donor_metadata.tsv")
readr::write_tsv(donor_counts_by_study, "revision_out/S9H_C3_donor_counts_by_study.tsv")
if (nrow(join_check) > 0) {
  readr::write_tsv(join_check, "revision_out/S9H_C3_annotation_join_check.tsv")
}

writeLines(
  c(
    "S9H (Lineage 3 cluster C3) donor-level trait analysis",
    "",
    "Consistency audit item: the analogous C2/C5 analysis was revised in",
    "response to R1-9. The same donor-level framework and multiple-testing",
    "correction are applied here so that both donor-trait tables in the",
    "manuscript are held to one standard.",
    "",
    "Donor/group identifiers:",
    "- GSE81608: source donor_id joined using Run = sample_name",
    "- GSE153855: source BioSample joined using Run = sample_name",
    "- GSE98887/GSE124742: metadata-defined fallback groups",
    "- Other datasets: existing orig.ident, prefixed by study",
    "",
    paste0("Total donor/groups analysed: ", dplyr::n_distinct(donor_meta$donor_id)),
    "",
    "Occupancy definition:",
    "number of a donor's Lineage 3 cells assigned to C3 divided by the total",
    "number of Lineage 3 cells contributed by that donor.",
    "",
    "Statistical tests:",
    "- BMI, age and HbA1c: Spearman correlation",
    "- Gender: Wilcoxon rank-sum",
    "- Ethnicity: Kruskal-Wallis",
    paste0("- Multiple testing: Benjamini-Hochberg across the ", nrow(results),
           " tests in this analysis (family = '", BH_FAMILY, "')"),
    "",
    "Supersedes the previous S9H table, which used raw orig.ident as the donor",
    "unit (collapsing each pooled GSE study into a single donor), parametric",
    "tests on bounded proportions, and no multiple-testing correction."
  ),
  "revision_out/S9H_C3_analysis_notes.txt"
)

capture.output(sessionInfo(), file = "revision_out/S9H_C3_sessionInfo.txt")

cat("\nDone. Outputs in revision_out/S9H_C3_*\n")
