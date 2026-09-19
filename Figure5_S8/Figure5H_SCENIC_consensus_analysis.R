################################################################################
# Project: scRNA beta cell atlas — Fig 5H SCENIC consensus analysis
# File: Fig5H_SCENIC_consensus_analysis.R
# Author: Elhadi Iich
#
# Description:
#   Performs donor-aware regulon testing for Lineage2 C3 using pySCENIC
#   consensus networks generated across 50 GRNBoost2 runs. The 70% consensus
#   network is used as the primary analysis, with 50% and 80% thresholds used
#   for sensitivity analysis.
################################################################################



library(dplyr)
library(tibble)
library(readr)
library(lme4)
library(SCopeLoomR)
library(emmeans)

CONSENSUS_PCT <- 70
GROUP_ORDER <- c("Lineage2_1", "Shared", "Lineage2_2")

################################################################################
### 1. Load SCENIC AUCell + metadata
################################################################################

loom_path <- sprintf(
  "pySCENIC/revision_out/scenic_consensus%d.loom",
  CONSENSUS_PCT
)

loom <- open_loom(loom_path)
regulonAUC_obj <- get_regulons_AUC(
  loom,
  column.attr.name = "RegulonsAUC"
)

regulonAUC_mat <- SummarizedExperiment::assay(
  regulonAUC_obj,
  "AUC"
)
close_loom(loom)

cellInfo <- readRDS("cellInfo_slingshotBranchID.Rds")

donor_obj <- readRDS(
  "./revision/Fig4_S6/lin2.beta.reintegrated.rds"
)

meta <- donor_obj@meta.data

common_cells <- Reduce(
  intersect,
  list(
    colnames(regulonAUC_mat),
    rownames(cellInfo),
    rownames(meta)
  )
)

cat("Matched cells:", length(common_cells), "/", ncol(regulonAUC_mat), "\n")

regulonAUC_mat <- regulonAUC_mat[, common_cells, drop = FALSE]

################################################################################
### 2. Load explicit donor mappings
################################################################################

ann81608 <- read_csv(
  "./revision/R1-9__R2-15/revision_out/GSE81608.csv",
  show_col_types = FALSE
) %>%
  transmute(
    run_id = as.character(Run),
    mapped_id = paste0("GSE81608::", donor_id)
  )

ann153855 <- read_csv(
  "./revision/R1-9__R2-15/revision_out/GSE153855.csv",
  show_col_types = FALSE
) %>%
  transmute(
    run_id = as.character(Run),
    mapped_id = paste0("GSE153855::", BioSample)
  )

################################################################################
### 3. Technical-run -> biological sample mappings
################################################################################

run_map <- tribble(
  ~run_id, ~mapped_id,
  
  ## GSE148073
  "SRR11477627", "GSE148073::SAMN14534753",
  "SRR11477628", "GSE148073::SAMN14534753",
  "SRR11477629", "GSE148073::SAMN14534752",
  "SRR11477630", "GSE148073::SAMN14534751",
  "SRR11477631", "GSE148073::SAMN14534750",
  "SRR11477632", "GSE148073::SAMN14534749",
  "SRR11477633", "GSE148073::SAMN14534748",
  "SRR11477634", "GSE148073::SAMN14534748",
  "SRR11477635", "GSE148073::SAMN14534747",
  "SRR11477636", "GSE148073::SAMN14534747",
  "SRR11477637", "GSE148073::SAMN14534747",
  "SRR11477638", "GSE148073::SAMN14534746",
  "SRR11477639", "GSE148073::SAMN14534746",
  "SRR11477640", "GSE148073::SAMN14534745",
  "SRR11477641", "GSE148073::SAMN14534745",
  "SRR11477642", "GSE148073::SAMN14534744",
  "SRR11477643", "GSE148073::SAMN14534744",
  "SRR11477644", "GSE148073::SAMN14534744",
  "SRR11477645", "GSE148073::SAMN14534744",
  "SRR11477646", "GSE148073::SAMN14534744",
  "SRR11477647", "GSE148073::SAMN14534744",
  "SRR11477648", "GSE148073::SAMN14534743",
  "SRR11477649", "GSE148073::SAMN14534743",
  "SRR11477650", "GSE148073::SAMN14534743",
  "SRR11477651", "GSE148073::SAMN14534742",
  "SRR11477652", "GSE148073::SAMN14534742",
  "SRR11477653", "GSE148073::SAMN14534741",
  "SRR11477654", "GSE148073::SAMN14534741",
  "SRR11477655", "GSE148073::SAMN14534741",
  "SRR11477656", "GSE148073::SAMN14534741",
  "SRR11477657", "GSE148073::SAMN14534740",
  "SRR11477658", "GSE148073::SAMN14534740",
  
  ## GSE150724
  "SRR11804772", "GSE150724::SAMN11633049",
  "SRR11804773", "GSE150724::SAMN11633049",
  "SRR11804774", "GSE150724::SAMN11963659",
  "SRR11804775", "GSE150724::SAMN11963659",
  "SRR11804776", "GSE150724::SAMN12227196",
  "SRR11804777", "GSE150724::SAMN12227196",
  
  ## GSE84133
  "SRR3879601", "GSE84133::SAMN05379702",
  "SRR3879602", "GSE84133::SAMN05379702",
  "SRR3879603", "GSE84133::SAMN05379702",
  "SRR3879604", "GSE84133::SAMN05379703",
  "SRR3879605", "GSE84133::SAMN05379703",
  "SRR3879606", "GSE84133::SAMN05379703",
  "SRR3879607", "GSE84133::SAMN05379704",
  "SRR3879608", "GSE84133::SAMN05379704",
  "SRR3879609", "GSE84133::SAMN05379704",
  "SRR3879610", "GSE84133::SAMN05379704",
  "SRR3879611", "GSE84133::SAMN05379705",
  "SRR3879612", "GSE84133::SAMN05379705",
  "SRR3879613", "GSE84133::SAMN05379705"
)

run_map <- bind_rows(
  run_map,
  ann81608,
  ann153855
) %>%
  distinct(run_id, .keep_all = TRUE)

################################################################################
### 4. Build donor IDs for the C3 cells
################################################################################

meta_c3 <- meta[common_cells, , drop = FALSE] %>%
  rownames_to_column("cell_id") %>%
  mutate(
    run_id = case_when(
      grepl("^SRR", orig.ident) ~ as.character(orig.ident),
      grepl("^SRR", sample_name) ~ as.character(sample_name),
      TRUE ~ as.character(sample_name)
    )
  ) %>%
  left_join(run_map, by = "run_id") %>%
  mutate(
    donor_id = case_when(
      
      Study %in% c(
        "GSE81608",
        "GSE153855",
        "GSE148073",
        "GSE150724",
        "GSE84133"
      ) ~ mapped_id,
      
      Study == "GSE114297" ~ paste0(
        "GSE114297::",
        orig.ident
      ),
      
      Study == "HPAP" ~ paste0(
        "HPAP::",
        orig.ident
      ),
      
      Study %in% c(
        "GSE98887",
        "GSE124742"
      ) ~ NA_character_,
      
      TRUE ~ NA_character_
    ),
    
    LineageGroup = factor(
      as.character(cellInfo[cell_id, "LineageGroup"]),
      levels = GROUP_ORDER
    )
  )

################################################################################
### 5. Check donor reconstruction
################################################################################

inferential_meta <- meta_c3 %>%
  filter(!Study %in% c("GSE98887", "GSE124742"))


cat("\nDonors represented by study:\n")

print(
  inferential_meta %>%
    dplyr::distinct(Study, donor_id) %>%
    dplyr::count(Study, name = "n_donors")
)

cat("\nExcluded from donor-level inference:\n")

print(
  meta_c3 %>%
    dplyr::filter(Study %in% c("GSE98887", "GSE124742")) %>%
    dplyr::count(Study, name = "n_cells")
)

################################################################################
### 6. Donor-aware regulon testing
################################################################################

regulon_names <- rownames(regulonAUC_mat)

test_regulon <- function(reg) {
  
  dat <- tibble(
    cell_id = common_cells,
    auc = as.numeric(regulonAUC_mat[reg, ])
  ) %>%
    left_join(
      meta_c3 %>%
        dplyr::select(cell_id, donor_id, LineageGroup),
      by = "cell_id"
    ) %>%
    dplyr::filter(
      !is.na(donor_id),
      !is.na(LineageGroup)
    ) %>%
    dplyr::group_by(
      donor_id,
      LineageGroup
    ) %>%
    dplyr::summarise(
      mean_auc = mean(auc, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (
    n_distinct(dat$donor_id) < 2 ||
    n_distinct(dat$LineageGroup) < 2
  ) {
    return(
      tibble(
        regulon = reg,
        chisq = NA_real_,
        p_value = NA_real_,
        peak_group = NA_character_,
        n_donors = n_distinct(dat$donor_id),
        singular = NA
      )
    )
  }
  
  full <- lme4::lmer(
    mean_auc ~ LineageGroup + (1 | donor_id),
    data = dat,
    REML = FALSE
  )
  
  singular_fit <- lme4::isSingular(full)
  
  null <- lme4::lmer(
    mean_auc ~ 1 + (1 | donor_id),
    data = dat,
    REML = FALSE
  )
  
  cmp <- anova(null, full)
  
  group_means <- tapply(
    dat$mean_auc,
    dat$LineageGroup,
    mean,
    na.rm = TRUE
  )
  
  tibble(
    regulon = reg,
    chisq = unname(cmp$Chisq[2]),
    p_value = cmp$`Pr(>Chisq)`[2],
    peak_group = names(group_means)[which.max(group_means)],
    n_donors = n_distinct(dat$donor_id),
    singular = singular_fit
  )
}

stats_tbl <- dplyr::bind_rows(
  lapply(regulon_names, test_regulon)
) %>%
  dplyr::filter(!is.na(p_value)) %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(p_adj)

################################################################################
### 7. Save results
################################################################################

readr::write_tsv(
  stats_tbl,
  "Fig5H_regulon_donor_mixed_model.tsv"
)

cat("\nTop results:\n")
print(head(stats_tbl, 20))

cat(
  "\nSignificant regulons, BH < 0.05:",
  sum(stats_tbl$p_adj < 0.05),
  "\n"
)

table(stats_tbl$singular)

stats_tbl %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::select(regulon, p_adj, peak_group, singular)


################################################################################
### 8. Sensitivity analysis across 50%, 70%, and 80% consensus thresholds
################################################################################

run_threshold_analysis <- function(pct) {
  
  cat("\n\n===== Consensus threshold:", pct, "% =====\n")
  
  loom_path_thr <- sprintf(
    "pySCENIC/revision_out/scenic_consensus%d.loom",
    pct
  )
  
  loom_thr <- open_loom(loom_path_thr)
  
  auc_obj_thr <- get_regulons_AUC(
    loom_thr,
    column.attr.name = "RegulonsAUC"
  )
  
  auc_mat_thr <- SummarizedExperiment::assay(
    auc_obj_thr,
    "AUC"
  )
  
  close_loom(loom_thr)
  
  cells_thr <- intersect(
    colnames(auc_mat_thr),
    meta_c3$cell_id
  )
  
  cat(
    "Matched cells:",
    length(cells_thr),
    "/",
    ncol(auc_mat_thr),
    "\n"
  )
  
  
  auc_mat_thr <- auc_mat_thr[, cells_thr, drop = FALSE]
  
  test_one <- function(reg) {
    
    dat <- tibble(
      cell_id = cells_thr,
      auc = as.numeric(auc_mat_thr[reg, ])
    ) %>%
      dplyr::left_join(
        meta_c3 %>%
          dplyr::select(
            cell_id,
            donor_id,
            LineageGroup
          ),
        by = "cell_id"
      ) %>%
      dplyr::filter(
        !is.na(donor_id),
        !is.na(LineageGroup)
      ) %>%
      dplyr::group_by(
        donor_id,
        LineageGroup
      ) %>%
      dplyr::summarise(
        mean_auc = mean(auc, na.rm = TRUE),
        .groups = "drop"
      )
    
    full <- lme4::lmer(
      mean_auc ~ LineageGroup + (1 | donor_id),
      data = dat,
      REML = FALSE
    )
    
    singular_fit <- lme4::isSingular(full)
    
    null <- lme4::lmer(
      mean_auc ~ 1 + (1 | donor_id),
      data = dat,
      REML = FALSE
    )
    
    cmp <- anova(null, full)
    
    group_means <- tapply(
      dat$mean_auc,
      dat$LineageGroup,
      mean,
      na.rm = TRUE
    )
    
    tibble(
      consensus_pct = pct,
      regulon = reg,
      chisq = unname(cmp$Chisq[2]),
      p_value = cmp$`Pr(>Chisq)`[2],
      peak_group = names(group_means)[which.max(group_means)],
      n_donors = dplyr::n_distinct(dat$donor_id),
      singular = singular_fit
    )
  }
  
  res <- dplyr::bind_rows(
    lapply(
      rownames(auc_mat_thr),
      test_one
    )
  ) %>%
    dplyr::mutate(
      p_adj = p.adjust(
        p_value,
        method = "BH"
      )
    ) %>%
    dplyr::arrange(p_adj)
  
  cat(
    "Regulons:",
    nrow(res),
    "\nSignificant at BH < 0.05:",
    sum(res$p_adj < 0.05),
    "\n"
  )
  
  return(res)
}


threshold_results <- dplyr::bind_rows(
  lapply(
    c(50, 70, 80),
    run_threshold_analysis
  )
)

readr::write_tsv(
  threshold_results,
  "Fig5H_regulon_threshold_sensitivity.tsv"
)


################################################################################
### 9. Compare significant 70% regulons across thresholds
################################################################################

sig70 <- stats_tbl %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::pull(regulon)

sig70_threshold_check <- threshold_results %>%
  dplyr::filter(regulon %in% sig70) %>%
  dplyr::select(
    consensus_pct,
    regulon,
    p_adj,
    peak_group,
    singular
  ) %>%
  dplyr::arrange(
    regulon,
    consensus_pct
  )

print(
  sig70_threshold_check,
  n = Inf
)

readr::write_tsv(
  sig70_threshold_check,
  "Fig5H_sig70_regulons_across_thresholds.tsv"
)

################################################################################
### 10. Summarise threshold robustness of 70% significant regulons
################################################################################

robustness_summary <- sig70_threshold_check %>%
  dplyr::group_by(regulon) %>%
  dplyr::summarise(
    n_thresholds_present = dplyr::n(),
    n_thresholds_FDR05 = sum(p_adj < 0.05),
    peak_groups = paste(unique(peak_group), collapse = ";"),
    same_peak_all = dplyr::n_distinct(peak_group) == 1,
    any_singular = any(singular),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_thresholds_FDR05),
    dplyr::desc(same_peak_all),
    regulon
  )

print(robustness_summary, n = Inf)

readr::write_tsv(
  robustness_summary,
  "Fig5H_regulon_robustness_summary.tsv"
)

################################################################################
### 11. Donor-level data for robust regulons in revised Fig 5H
################################################################################

robust_regulons <- robustness_summary %>%
  dplyr::filter(
    n_thresholds_FDR05 >= 2,
    same_peak_all,
    !any_singular
  ) %>%
  dplyr::pull(regulon)

cat("\nRobust regulons retained for revised Fig 5H:\n")
print(robust_regulons)

fig5h_donor_data <- lapply(
  robust_regulons,
  function(reg) {
    
    tibble(
      cell_id = common_cells,
      auc = as.numeric(regulonAUC_mat[reg, ])
    ) %>%
      dplyr::left_join(
        meta_c3 %>%
          dplyr::select(
            cell_id,
            Study,
            donor_id,
            LineageGroup
          ),
        by = "cell_id"
      ) %>%
      dplyr::filter(
        !is.na(donor_id),
        !is.na(LineageGroup)
      ) %>%
      dplyr::group_by(
        Study,
        donor_id,
        LineageGroup
      ) %>%
      dplyr::summarise(
        mean_auc = mean(auc, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        regulon = reg
      )
  }
) %>%
  dplyr::bind_rows()

readr::write_tsv(
  fig5h_donor_data,
  "Fig5H_regulon_donor_data.tsv"
)

print(
  fig5h_donor_data %>%
    dplyr::group_by(regulon, LineageGroup) %>%
    dplyr::summarise(
      n_donors = dplyr::n_distinct(donor_id),
      mean_auc = mean(mean_auc),
      .groups = "drop"
    ),
  n = Inf
)

################################################################################
### 12. Pairwise lineage comparisons for robust Fig 5H regulons
################################################################################


pairwise_results <- lapply(
  robust_regulons,
  function(reg) {
    
    dat <- fig5h_donor_data %>%
      dplyr::filter(regulon == reg)
    
    fit <- lme4::lmer(
      mean_auc ~ LineageGroup + (1 | donor_id),
      data = dat,
      REML = FALSE
    )
    
    emm <- emmeans::emmeans(
      fit,
      ~ LineageGroup
    )
    
    contrasts <- as.data.frame(
      emmeans::contrast(
        emm,
        method = "pairwise",
        adjust = "none"
      )
    )
    
    contrasts %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        regulon = reg
      ) %>%
      dplyr::select(
        regulon,
        contrast,
        estimate,
        SE,
        df,
        t.ratio,
        p.value
      )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_adj = p.adjust(
      p.value,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(
    regulon,
    p_adj
  )

print(
  pairwise_results,
  n = Inf
)

readr::write_tsv(
  pairwise_results,
  "Fig5H_pairwise_mixed_model.tsv"
)

################################################################################
### 13. Supplementary Table S16
################################################################################

# Primary 70% results: retain the 19 significant regulons
table_s16 <- stats_tbl %>%
  dplyr::filter(p_adj < 0.05) %>%
  dplyr::transmute(
    regulon,
    chisq_70 = chisq,
    p_value_70 = p_value,
    p_adj_70 = p_adj,
    peak_group_70 = peak_group,
    n_donors,
    singular_70 = singular
  )

threshold_wide <- sig70_threshold_check %>%
  dplyr::filter(consensus_pct %in% c(50, 80)) %>%
  tidyr::pivot_wider(
    names_from = consensus_pct,
    values_from = c(p_adj, peak_group, singular),
    names_glue = "{.value}_{consensus_pct}"
  )

robust_tbl <- robustness_summary %>%
  dplyr::mutate(
    robust = n_thresholds_FDR05 >= 2 &
      same_peak_all &
      !any_singular
  )

table_s16 <- table_s16 %>%
  dplyr::left_join(threshold_wide, by = "regulon") %>%
  dplyr::left_join(robust_tbl, by = "regulon") %>%
  dplyr::arrange(
    dplyr::desc(robust),
    p_adj_70
  )

readr::write_tsv(
  table_s16,
  "table_s16_SCENIC_regulon_robustness.tsv"
)

print(table_s16, n = Inf)

