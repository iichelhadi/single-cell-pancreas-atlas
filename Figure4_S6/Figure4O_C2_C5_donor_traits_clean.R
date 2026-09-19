################################################################################
# Figure 4O — donor-level C2/C5 trait associations
# B9 / C15
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
})

dir.create("revision_out", showWarnings = FALSE)

################################################################################
# Load Lineage 2 object
################################################################################

lin2_obj <- readRDS(
  "lin2.beta.reintegrated.rds"
)


################################################################################
# Extract current cell metadata
################################################################################

metadata <- lin2_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::transmute(
    cell_id,
    orig.ident = as.character(orig.ident),
    sample_name = as.character(sample_name),
    Study = as.character(Study),
    cluster = as.character(monocle_clusters2),
    Age = as.numeric(Age),
    BMI = as.numeric(BMI),
    Hba1c = as.numeric(Hba1c),
    Gender = as.character(Gender),
    Ethnicity = as.character(Ethnicity),
    Condition = as.character(Condition)
  )


################################################################################
# Read source annotation mappings
################################################################################

ann81608 <- readr::read_csv(
  "./revision_out/GSE81608.csv",
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    sample_name = as.character(Run),
    annotated_donor = as.character(donor_id)
  ) %>%
  dplyr::filter(
    !is.na(sample_name),
    sample_name != "",
    !is.na(annotated_donor),
    annotated_donor != ""
  ) %>%
  dplyr::distinct()

ann153855 <- readr::read_csv(
  "./revision_out/GSE153855.csv",
  show_col_types = FALSE
) %>%
  dplyr::transmute(
    sample_name = as.character(Run),
    annotated_biosample = as.character(BioSample)
  ) %>%
  dplyr::filter(
    !is.na(sample_name),
    sample_name != "",
    !is.na(annotated_biosample),
    annotated_biosample != ""
  ) %>%
  dplyr::distinct()

################################################################################
# Verify that each Run has only one source identifier
################################################################################

conflict81608 <- ann81608 %>%
  dplyr::count(sample_name, name = "n_mappings") %>%
  dplyr::filter(n_mappings > 1)

conflict153855 <- ann153855 %>%
  dplyr::count(sample_name, name = "n_mappings") %>%
  dplyr::filter(n_mappings > 1)


################################################################################
# Join source identifiers to current cells
################################################################################

metadata <- metadata %>%
  dplyr::left_join(
    ann81608,
    by = "sample_name"
  ) %>%
  dplyr::left_join(
    ann153855,
    by = "sample_name"
  )

################################################################################
# Verify annotation joins
################################################################################

join_check <- metadata %>%
  dplyr::filter(
    orig.ident %in% c("GSE81608", "GSE153855")
  ) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    
    n_matched = dplyr::case_when(
      dplyr::first(orig.ident) == "GSE81608" ~
        sum(!is.na(annotated_donor)),
      
      dplyr::first(orig.ident) == "GSE153855" ~
        sum(!is.na(annotated_biosample))
    ),
    
    n_identifiers = dplyr::case_when(
      dplyr::first(orig.ident) == "GSE81608" ~
        dplyr::n_distinct(annotated_donor, na.rm = TRUE),
      
      dplyr::first(orig.ident) == "GSE153855" ~
        dplyr::n_distinct(annotated_biosample, na.rm = TRUE)
    ),
    
    .groups = "drop"
  )



################################################################################
# Construct fallback metadata signatures
################################################################################

fallback_studies <- c(
  "GSE98887",
  "GSE124742"
)

signature_value <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x
}

metadata <- metadata %>%
  dplyr::mutate(
    age_id = signature_value(Age),
    bmi_id = signature_value(BMI),
    hba1c_id = signature_value(Hba1c),
    gender_id = signature_value(Gender),
    ethnicity_id = signature_value(Ethnicity),
    condition_id = signature_value(Condition),
    
    donor_signature = paste(
      age_id,
      bmi_id,
      hba1c_id,
      gender_id,
      ethnicity_id,
      condition_id,
      sep = "|"
    )
  )

################################################################################
# Create final harmonized donor ID
################################################################################

metadata <- metadata %>%
  dplyr::mutate(
    donor_id = dplyr::case_when(
      
            orig.ident == "GSE81608" ~ paste0(
        "GSE81608::",
        annotated_donor
      ),
      
            orig.ident == "GSE153855" ~ paste0(
        "GSE153855::",
        annotated_biosample
      ),
      
            orig.ident %in% fallback_studies ~ paste0(
        orig.ident,
        "::metadata::",
        donor_signature
      ),
      
            TRUE ~ paste0(
        Study,
        "::",
        orig.ident
      )
    ),
    
    donor_id_source = dplyr::case_when(
      orig.ident == "GSE81608" ~
        "annotated donor ID",
      
      orig.ident == "GSE153855" ~
        "annotated BioSample",
      
      orig.ident %in% fallback_studies ~
        "metadata-defined group",
      
      TRUE ~
        "existing donor/sample identifier"
    )
  )


################################################################################
# Inspect final donor counts
################################################################################

donor_counts_by_study <- metadata %>%
  dplyr::group_by(
    Study,
    orig.ident,
    donor_id_source
  ) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    n_donor_groups = dplyr::n_distinct(donor_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Study, orig.ident)


################################################################################
# Verify metadata consistency within final donor IDs
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
    n_studies > 1 |
      n_age > 1 |
      n_BMI > 1 |
      n_HbA1c > 1 |
      n_Gender > 1 |
      n_Ethnicity > 1 |
      n_Condition > 1
  )



################################################################################
# Create one metadata row per donor/group
################################################################################

first_numeric <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  x[[1]]
}

first_character <- function(x) {
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
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



################################################################################
# Calculate donor-level C2 and C5 occupancy
################################################################################

donor_totals <- metadata %>%
  dplyr::count(
    donor_id,
    name = "total_lineage2_cells"
  )

donor_props <- metadata %>%
  dplyr::filter(cluster %in% c("2", "5")) %>%
  dplyr::count(
    donor_id,
    cluster,
    name = "n_cells"
  ) %>%
  tidyr::complete(
    donor_id = unique(metadata$donor_id),
    cluster = c("2", "5"),
    fill = list(n_cells = 0)
  ) %>%
  dplyr::left_join(
    donor_totals,
    by = "donor_id"
  ) %>%
  dplyr::mutate(
    proportion = n_cells / total_lineage2_cells
  ) %>%
  dplyr::left_join(
    donor_meta,
    by = "donor_id"
  ) %>%
  dplyr::arrange(
    donor_id,
    cluster
  )



################################################################################
# Statistical test functions
################################################################################

continuous_test <- function(dat, variable) {
  
  d <- dat %>%
    dplyr::filter(
      !is.na(.data[[variable]]),
      !is.na(proportion)
    )
  
  fit <- stats::cor.test(
    d[[variable]],
    d$proportion,
    method = "spearman",
    exact = FALSE
  )
  
  tibble::tibble(
    trait = variable,
    test_method = "Spearman correlation",
    n_donors = nrow(d),
    statistic = unname(fit$estimate),
    statistic_name = "rho",
    p_value = fit$p.value,
    group_summary = NA_character_
  )
}

sex_test <- function(dat) {
  
  d <- dat %>%
    dplyr::filter(
      Gender %in% c("male", "female"),
      !is.na(proportion)
    )
  
  
  fit <- stats::wilcox.test(
    proportion ~ Gender,
    data = d,
    exact = FALSE
  )
  
  summary_by_group <- d %>%
    dplyr::group_by(Gender) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_proportion = median(
        proportion,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  tibble::tibble(
    trait = "Gender",
    test_method = "Wilcoxon rank-sum",
    n_donors = dplyr::n_distinct(d$donor_id),
    statistic = unname(fit$statistic),
    statistic_name = "W",
    p_value = fit$p.value,
    group_summary = paste0(
      summary_by_group$Gender,
      ": n=",
      summary_by_group$n,
      ", median=",
      round(summary_by_group$median_proportion, 4),
      collapse = "; "
    )
  )
}

ethnicity_test <- function(dat) {
  
  d <- dat %>%
    dplyr::filter(
      !is.na(Ethnicity),
      Ethnicity != "",
      !is.na(proportion)
    )
  
  
  fit <- stats::kruskal.test(
    proportion ~ Ethnicity,
    data = d
  )
  
  tibble::tibble(
    trait = "Ethnicity",
    test_method = "Kruskal-Wallis",
    n_donors = dplyr::n_distinct(d$donor_id),
    statistic = unname(fit$statistic),
    statistic_name = "chi-square",
    p_value = fit$p.value,
    group_summary = NA_character_
  )
}

################################################################################
# Run C2 and C5 association analyses
################################################################################

results <- purrr::map_dfr(
  c("2", "5"),
  function(cl) {
    
    d <- donor_props %>%
      dplyr::filter(cluster == cl)
    
    dplyr::bind_rows(
      continuous_test(d, "BMI"),
      continuous_test(d, "Age"),
      continuous_test(d, "Hba1c"),
      sex_test(d),
      ethnicity_test(d)
    ) %>%
      dplyr::mutate(
        cluster = paste0("C", cl),
        .before = 1
      )
  }
) %>%
  dplyr::mutate(
    p_adjusted_BH = stats::p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(
    cluster,
    match(
      trait,
      c(
        "BMI",
        "Age",
        "Hba1c",
        "Gender",
        "Ethnicity"
      )
    )
  )

print(results, n = Inf)

################################################################################
# Save outputs
################################################################################

readr::write_tsv(
  results,
  "revision_out/C2_C5_donor_trait_associations.tsv"
)

readr::write_tsv(
  donor_props,
  "revision_out/C2_C5_donor_level_occupancy.tsv"
)

readr::write_tsv(
  donor_meta,
  "revision_out/C2_C5_donor_metadata.tsv"
)

readr::write_tsv(
  donor_counts_by_study,
  "revision_out/C2_C5_donor_counts_by_study.tsv"
)

readr::write_tsv(
  join_check,
  "revision_out/C2_C5_annotation_join_check.tsv"
)


################################################################################
### Revised Figure 4O table from corrected donor-level analysis
################################################################################

stats_table <- results %>%
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
      ),
      TRUE ~ sprintf("%.3f", statistic)
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

stats_table


library(flextable)
library(officer)

fig4O_table <- stats_table %>%
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
    i = seq(1, nrow(stats_table), by = 2),
    bg = "#FFFFFF",
    part = "body"
  ) %>%
  flextable::bg(
    i = seq(2, nrow(stats_table), by = 2),
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
fig4O_table
