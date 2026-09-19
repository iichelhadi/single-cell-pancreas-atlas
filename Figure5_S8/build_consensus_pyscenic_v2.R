###############################################################################
# Consensus GRN builder for pySCENIC (reviewer C9) — v2
# CORRECTED: reads/writes WITH headers (TF, target, importance), matching
# pyscenic's own adjacency file convention (see run_grnboost2_multi_v2.py
# header comment for why headerless was wrong). Also adds an actual
# threshold-comparison step so 50/70/80% is a real sensitivity analysis,
# not three unexamined output files.
#
# INPUTS:
#   grn_runs/grn_output_run*.tsv   N adjacency lists, WITH header
#                                  (TF, target, importance)
#
# OUTPUTS (revision_out/):
#   grn_edge_recurrence.tsv        every edge's frequency + mean/median weight
#   consensus_grn_50.tsv           TF  target  importance   (WITH header)
#   consensus_grn_70.tsv
#   consensus_grn_80.tsv
#   threshold_comparison.tsv       edge counts per threshold (regulon/target
#                                  overlap filled in after ctx/aucell — see
#                                  compare_thresholds_post_ctx.R)
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})

dir.create("revision_out", showWarnings = FALSE)

THRESHOLDS <- c(0.50, 0.70, 0.80)

run_files <- list.files("grn_runs", pattern = "^grn_output_run.*\\.tsv$",
                         full.names = TRUE)

cat("Aggregating", length(run_files), "GRNBoost2 runs\n")

read_run <- function(f) {
  d <- fread(f, header = TRUE)   # v2: files now have real headers
  d[, .(TF, target, importance)]
}
all_runs <- rbindlist(lapply(run_files, read_run), idcol = "run")
n_runs   <- length(run_files)

edge_stats <- all_runs[, .(
  n_runs_present = uniqueN(run),
  freq           = uniqueN(run) / n_runs,
  mean_importance   = mean(importance),
  median_importance = median(importance)
), by = .(TF, target)]
fwrite(edge_stats, "revision_out/grn_edge_recurrence.tsv", sep = "\t")

cat(sprintf("Total unique edges observed: %d\n", nrow(edge_stats)))

threshold_summary <- data.table(threshold_pct = numeric(), n_edges = integer(),
                                 n_unique_TFs = integer(), n_unique_targets = integer())

for (th in THRESHOLDS) {
  consensus <- edge_stats[freq >= th]
  pct <- round(100 * th)
  out_file <- sprintf("revision_out/consensus_grn_%d.tsv", pct)

  out_df <- consensus[, .(TF, target, importance = mean_importance)]
  fwrite(out_df, out_file, sep = "\t", col.names = TRUE, quote = FALSE)

  threshold_summary <- rbind(threshold_summary, data.table(
    threshold_pct = pct,
    n_edges = nrow(consensus),
    n_unique_TFs = uniqueN(consensus$TF),
    n_unique_targets = uniqueN(consensus$target)
  ))

  cat(sprintf("Threshold >=%d%%: %d / %d edges retained (%d TFs, %d targets) -> %s\n",
              pct, nrow(consensus), nrow(edge_stats),
              uniqueN(consensus$TF), uniqueN(consensus$target), out_file))
}

## --- pairwise Jaccard overlap of edge sets between thresholds --------------
edge_sets <- lapply(THRESHOLDS, function(th) {
  edge_stats[freq >= th, paste(TF, target, sep = "|")]
})
names(edge_sets) <- paste0(round(100 * THRESHOLDS), "%")

## --- split-half stability: first 25 runs vs second 25 runs -----------------

run_ids <- sort(unique(all_runs$run))


half1_ids <- run_ids[1:25]
half2_ids <- run_ids[26:50]

build_half_stats <- function(run_subset) {
  x <- all_runs[run %in% run_subset]

  x[, .(
    n_runs_present = uniqueN(run),
    freq = uniqueN(run) / length(run_subset),
    mean_importance = mean(importance),
    median_importance = median(importance)
  ), by = .(TF, target)]
}

half1_stats <- build_half_stats(half1_ids)
half2_stats <- build_half_stats(half2_ids)

split_summary <- data.table()

for (th in THRESHOLDS) {
  pct <- round(100 * th)

  set1 <- half1_stats[
    freq >= th,
    paste(TF, target, sep = "|")
  ]

  set2 <- half2_stats[
    freq >= th,
    paste(TF, target, sep = "|")
  ]

  union_edges <- union(set1, set2)

  split_jaccard <- if (length(union_edges) == 0) {
    NA_real_
  } else {
    length(intersect(set1, set2)) / length(union_edges)
  }

  split_summary <- rbind(
    split_summary,
    data.table(
      threshold_pct = pct,
      half1_edges = length(set1),
      half2_edges = length(set2),
      shared_edges = length(intersect(set1, set2)),
      union_edges = length(union_edges),
      split_half_jaccard = split_jaccard
    )
  )
}

fwrite(
  split_summary,
  "revision_out/split_half_stability.tsv",
  sep = "\t"
)

cat("\nSplit-half stability:\n")
print(split_summary)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
pairs <- combn(names(edge_sets), 2, simplify = FALSE)
jaccard_tbl <- rbindlist(lapply(pairs, function(p) {
  data.table(threshold_A = p[1], threshold_B = p[2],
             jaccard_edge_overlap = jaccard(edge_sets[[p[1]]], edge_sets[[p[2]]]))
}))

fwrite(threshold_summary, "revision_out/threshold_comparison.tsv", sep = "\t")
fwrite(jaccard_tbl, "revision_out/threshold_jaccard_overlap.tsv", sep = "\t")

cat("\nEdge-count summary:\n"); print(threshold_summary)
cat("\nJaccard edge overlap between thresholds:\n"); print(jaccard_tbl)
cat("\nNOTE: this covers edge-level overlap only. After running ctx/aucell for\n")
cat("each threshold, also compare regulon-target-set overlap and AUCell score\n")
cat("correlation for shared regulons (GPT's point #4) — do this once .loom\n")
cat("files exist for at least two thresholds; see compare_thresholds_post_ctx.R.\n")

## ===========================================================================
## Environment snapshot (R side). Conda/pyscenic side is captured separately
## in run_pyscenic_ctx_aucell.pbs, in the same job as the ctx/aucell run.
## ===========================================================================
writeLines(capture.output(sessionInfo()), "revision_out/sessionInfo_consensus.txt")
ns <- loadedNamespaces()
ver <- data.frame(package = ns,
                   version = vapply(ns, function(p) as.character(utils::packageVersion(p)), ""))
write.table(ver[order(ver$package), ], "revision_out/package_versions_consensus.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nBEFORE the full ctx/aucell run: test on a small subset.\n")
cat("  head -5 revision_out/consensus_grn_70.tsv > revision_out/test_subset.tsv\n")
cat("  pyscenic ctx revision_out/test_subset.tsv <feather_dbs> \\\n")
cat("    --expression_mtx_fname <small_expr_subset.csv> --output test_out.csv\n")
cat("A column-format mismatch produces silently empty/zero regulons, not an\n")
cat("error — cheap to test now, expensive to discover after the full walltime.\n")
