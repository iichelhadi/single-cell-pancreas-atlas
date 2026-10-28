#!/usr/bin/env Rscript
# Build the merged scRNA-seq atlas: per-sample QC + DoubletFinder, merge,
# attach donor metadata, drop excluded studies.

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ---- CONFIG -----------------------------------------------------------------
data_root     <- "/path/to/data"                 # parent of the scRNAseq/ folder
samples_file  <- "samples.tsv"                   # one row per sample (see README)
donor_meta    <- "meta_data.txt"        # donor metadata, needs sample_id column
out_dir       <- "output"
exclude_studies <- c("GSE162708",
                     "GSE167880",
                     "PRJNA594252",
                     "GSE110154",
                     "GSE150212")
qc_vars       <- c("percent.mt",
                   "log10GenesPerUMI",
                   "nFeature_RNA",
                   "nCount_RNA")
# -----------------------------------------------------------------------------

set.seed(123)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
samples <- read.delim(samples_file,
                      check.names = FALSE,
                      stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(samples$sample_id),
          !anyDuplicated(samples$merge_prefix))

# Protein-coding genes only
gene_types <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"),
  columns = c("SYMBOL", "GENETYPE"),
  keytype = "SYMBOL"
)
coding_symbols <- unique(gene_types$SYMBOL[gene_types$GENETYPE == "protein-coding"])

resolve <- function(p)
  if (grepl("^/", p))
    p
else
  file.path(data_root, p)

process_sample <- function(s) {
  message("Processing ", s$sample_id, " (", s$merge_prefix, ")")
  
  files <- vapply(s[c("matrix", "features", "barcodes")], resolve, character(1))
  missing <- files[!file.exists(files)]
  if (length(missing))
    stop("Missing input(s) for ",
         s$sample_id,
         ": ",
         paste(missing, collapse = ", "))
  
  counts <- ReadMtx(mtx = files[["matrix"]],
                    features = files[["features"]],
                    cells = files[["barcodes"]])
  obj <- CreateSeuratObject(
    counts,
    min.cells = 0,
    min.features = as.integer(s$min_cells_genes),
    project = s$project
  )
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj <- obj[rownames(obj) %in% coding_symbols, ]
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  
  # Per-sample QC filter from the manifest; restricted to known QC columns.
  expr <- parse(text = s$qc_expression)[[1]]
  bad <- setdiff(all.names(expr, functions = FALSE, unique = TRUE), qc_vars)
  if (length(bad))
    stop("Unexpected QC variable for ",
         s$sample_id,
         ": ",
         paste(bad, collapse = ","))
  keep <- eval(expr, envir = as.list(obj@meta.data), enclos = baseenv())
  if (length(keep) != ncol(obj) ||
      anyNA(keep))
    stop("Invalid QC mask for ", s$sample_id)
  obj <- subset(obj, cells = colnames(obj)[keep])
  if (ncol(obj) < 3L)
    stop("Too few cells after QC: ", s$sample_id)
  
  # Normalized copy used only to call doublets; the un-normalized object is merged.
  raw <- obj
  dims <- seq_len(as.integer(s$pcs))
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  hvg <- VariableFeatures(obj)
  obj <- ScaleData(
    obj,
    features = hvg,
    vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt")
  )
  obj <- RunPCA(obj, features = hvg, npcs = 50)
  obj <- FindNeighbors(obj, dims = dims)
  obj <- FindClusters(obj, resolution = 0.1)
  obj <- RunUMAP(obj, dims = dims)
  
  sweep <- paramSweep_v3(obj, PCs = dims, sct = FALSE)
  bcmvn <- find.pK(summarizeSweep(sweep, GT = FALSE))
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic <- modelHomotypic(obj$seurat_clusters)
  nExp <- round(round(as.numeric(s$doublet_rate) * ncol(obj)) * (1 - homotypic))
  obj <- doubletFinder_v3(
    obj,
    PCs = dims,
    pN = 0.25,
    pK = pK,
    nExp = nExp,
    reuse.pANN = FALSE,
    sct = FALSE
  )
  
  # Drop the pANN column (second to last) and rename the classification column.
  obj@meta.data[ncol(obj@meta.data) - 1L] <- NULL
  names(obj@meta.data)[ncol(obj@meta.data)] <- "doublet_status"
  raw@meta.data <- obj@meta.data
  raw
}

processed <- lapply(seq_len(nrow(samples)), function(i)
  process_sample(samples[i, ]))
merged <- merge(processed[[1]], processed[-1], add.cell.ids = samples$merge_prefix)
rm(processed)
gc()
saveRDS(merged, file.path(out_dir, "merged_seurat_raw.rds"))

# study = first barcode field, sample_name = second
barcodes <- rownames(merged@meta.data)
parts <- strsplit(barcodes, "_", fixed = TRUE)
merged$study       <- vapply(parts, `[`, character(1), 1L)
merged$sample_name <- vapply(parts, function(x)
  if (length(x) >= 2L)
    x[2]
  else
    "", character(1))

# Join donor metadata by sample_name
meta <- read.delim(donor_meta,
                   check.names = FALSE,
                   stringsAsFactors = FALSE)
if (!"sample_id" %in% names(meta))
  stop("Donor metadata requires a sample_id column")
if (anyDuplicated(meta$sample_id))
  stop("Duplicate sample_id in donor metadata")
extra <- meta[match(merged$sample_name, meta$sample_id), setdiff(names(meta), "sample_id"), drop = FALSE]
rownames(extra) <- barcodes
merged@meta.data <- cbind(merged@meta.data, extra)

write.table(
  merged@meta.data,
  file.path(out_dir, "meta_raw.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(merged, file.path(out_dir, "merged_seurat_raw2.rds"))

# Remove excluded studies (not a singlet filter)
tissue <- subset(merged, subset = !(study %in% exclude_studies))
write.table(
  tissue@meta.data,
  file.path(out_dir, "tissue_meta.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sce <- scDblFinder(
  as.SingleCellExperiment(tissue),
  samples = "orig.ident",
  BPPARAM = MulticoreParam(2),
  clusters = TRUE
)

tissue@meta.data <- cbind(tissue@meta.data, scDblFinder.class = sce@colData[, "scDblFinder.class"])
tissue@meta.data$doublet_final <- ifelse(
  tissue@meta.data$scDblFinder.class == "doublet" |
    tissue@meta.data$doublet_status == "Doublet",
  "Doublet",
  "Singlet"
)

saveRDS(tissue, file.path(out_dir, "scRNA_consolidation_part1.rds"))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
