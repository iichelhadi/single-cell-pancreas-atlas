# scRNA-seq object construction (CODE-1)

This is a manifest-driven refactor of `scRNA_consolidation_part1v20.R`, **not a new analytical workflow**. The accompanying `samples.tsv` was extracted from the source and retains its **116 sample entries in the original merge order**, their exact barcode prefixes, per-sample QC expressions, PCA dimensions, and expected doublet fractions. The 116 entries are the sample objects in this particular source script, not a claim about the total number of datasets in the full atlas.

## Run

```bash
export SCRNA_DATA_ROOT=/path/to/data-root
Rscript build_atlas.R samples.tsv /path/to/20231116_meta_data.txt /path/to/output
```

`SCRNA_DATA_ROOT` is the directory containing `scRNAseq/`. Input paths in the manifest are relative to this root. Source data and donor metadata are not included.

## Behavior retained

- Read original SoupX-corrected or platform-specific count matrices, as specified per row.
- Create Seurat objects, calculate mitochondrial/ribosomal percentages, restrict to protein-coding symbols, and apply **sample-specific** QC filters.
- Run normalization, 2,000 variable features, regression, 50-PC PCA, sample-specific neighbor/UMAP dimensions, and DoubletFinder parameter sweep.
- Preserve the source-specific expected doublet fractions (including non-8% values).
- Transfer doublet classifications to the **pre-normalization** Seurat objects; merge those objects in original order with original prefixes.
- Add metadata and save the three original RDS outputs and metadata TSV files.

## Validation limits

This refactor has been statically checked against the source, **not executed on the original count matrices**. Identical numerical output is not guaranteed without testing using the original R/Seurat/DoubletFinder versions, package databases, random seeds, and donor metadata. In particular, UMAP and clustering can change across versions. The original script uses `min.genes = 200`; this refactor maps it to `min.features = 200`, the corresponding documented Seurat argument. The source `doubletFinder_v3` and `paramSweep_v3` APIs may require the historical DoubletFinder version.

The source does **not** remove doublets at this stage: it retains cells and stores `doublet_status`. No `doublet_final` column is created here. That downstream decision remains outside the scope of this script.

Before claiming exact equivalence, compare cell barcodes, gene order, per-sample QC pass counts, doublet labels, and merged metadata against the original outputs. The manifest is source-derived; inspect it for the intended input files on your system.
