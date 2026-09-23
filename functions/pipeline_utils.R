# ============================================================================
# pipeline_utils.R
# ============================================================================
#
# Reusable helper functions for bulk RNA-seq analysis of Plasmidsaurus data.
#
# Source this file from any analysis notebook with:
#   source("pipeline_utils.R")
#
# Functions provided (in order of typical use):
#
#   --- Loading ---
#   load_plasmidsaurus()        - Read and parse Plasmidsaurus expression matrix
#
#   --- QC (run before analysis) ---
#   load_mapping_stats()        - Read Plasmidsaurus mapping stats CSV
#   compute_qc_metrics()        - Per-sample library size, detection, biotype %, etc.
#   check_sex_markers()         - Verify sex assignment using sex-linked gene expression
#   plot_qc_summary()           - Multi-panel QC figure
#   plot_sample_correlation()   - Pearson correlation heatmap between samples
#   flag_qc_outliers()          - Apply thresholds and flag samples for review
#   run_qc_report()             - One-call wrapper that does all of the above
#
#   --- Differential expression ---
#   build_dds()                 - Build and fit a DESeq2 object
#   extract_contrast()          - Get a clean results table for one contrast
#   summarize_significance()    - Print significant gene counts at standard thresholds
#
#   --- Plotting ---
#   plot_pca_custom()           - PCA colored by any metadata column(s)
#   plot_sample_distance()      - Sample-to-sample distance heatmap
#   plot_volcano_custom()       - Volcano plot with sensible defaults
#
#   --- Output ---
#   write_results_excel()       - Save results table as Excel with formatting
#   write_outputs_for_enrichment() - Save full + significant + ranked outputs
#
# Dependencies (all should be in the bulkrnaseq conda environment):
#   DESeq2, tidyverse, pheatmap, RColorBrewer, ggrepel, openxlsx,
#   apeglm, ashr, EnhancedVolcano
#
# ============================================================================


# --------------------------------------------------------------------------
# load_plasmidsaurus()
# --------------------------------------------------------------------------
#' Load a Plasmidsaurus expression matrix
#'
#' Plasmidsaurus delivers a TSV with both CPM (normalized) and count (raw)
#' columns interleaved per sample, plus gene_id, gene_name, gene_biotype.
#' This function pulls out the count columns, rounds to integers (Plasmidsaurus
#' uses fractional multi-mapping assignment), and returns a clean counts matrix
#' plus a separate gene annotation table.
#'
#' @param tsv_path Path to the Plasmidsaurus expression matrix TSV file.
#' @param strip_prefix Optional regex pattern to strip from sample names
#'   (e.g., "^SS8VZG_" to remove the order prefix). Default NULL keeps names as-is.
#'
#' @return A list with three elements:
#'   - counts: integer matrix (genes x samples), rownames = Ensembl gene IDs
#'   - gene_info: data.frame with gene_id, gene_name, gene_biotype
#'   - sample_ids: character vector of sample IDs as they appear in the file
load_plasmidsaurus <- function(tsv_path, strip_prefix = NULL) {

  if (!file.exists(tsv_path)) {
    stop("File not found: ", tsv_path)
  }

  cat("Reading", tsv_path, "...\n")
  raw <- read.delim(tsv_path, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE)

  # Validate expected structure
  required_cols <- c("gene_id", "gene_name", "gene_biotype")
  missing <- setdiff(required_cols, colnames(raw))
  if (length(missing) > 0) {
    stop("Expected columns missing from file: ", paste(missing, collapse = ", "))
  }

  # Separate gene annotation
  gene_info <- raw[, required_cols]

  # Identify count columns (end in "_count") and CPM columns (end in "_cpm")
  count_cols <- grep("_count$", colnames(raw), value = TRUE)
  cpm_cols   <- grep("_cpm$",   colnames(raw), value = TRUE)

  if (length(count_cols) == 0) {
    stop("No '_count' columns found. Is this a Plasmidsaurus expression matrix?")
  }

  # Build counts matrix
  counts <- as.matrix(raw[, count_cols])
  rownames(counts) <- raw$gene_id

  # Strip "_count" suffix from sample names
  colnames(counts) <- sub("_count$", "", colnames(counts))

  # Optionally strip a prefix (e.g., the order ID like "SS8VZG_")
  if (!is.null(strip_prefix)) {
    colnames(counts) <- sub(strip_prefix, "", colnames(counts))
  }

  # Round fractional counts to integers (DESeq2 requires integers)
  n_fractional <- sum(counts != floor(counts))
  if (n_fractional > 0) {
    cat("  Rounding", n_fractional, "fractional count values to integers\n")
    counts <- round(counts)
  }
  storage.mode(counts) <- "integer"

  # Report
  cat("  Genes:  ", nrow(counts), "\n")
  cat("  Samples:", ncol(counts), "\n")
  cat("  Total CPM columns ignored:", length(cpm_cols), "\n")
  cat("  Library size range: ",
      format(min(colSums(counts)), big.mark = ","), "to",
      format(max(colSums(counts)), big.mark = ","), "reads\n")

  list(
    counts     = counts,
    gene_info  = gene_info,
    sample_ids = colnames(counts)
  )
}


# --------------------------------------------------------------------------
# build_dds()
# --------------------------------------------------------------------------
#' Build and fit a DESeq2 object
#'
#' Wraps DESeqDataSetFromMatrix + filtering + DESeq() into a single call.
#' Filters out genes with fewer than `min_count` reads in at least
#' `smallest_group_size` samples.
#'
#' @param counts Integer counts matrix (genes x samples).
#' @param metadata data.frame with one row per sample. rownames must match
#'   colnames(counts) exactly.
#' @param design A formula, e.g. ~ Exposure or ~ Sex + Exposure + Sex:Exposure.
#' @param min_count Minimum count for a gene to "pass" in a single sample. Default 10.
#' @param smallest_group_size Minimum number of samples that must pass for a gene
#'   to be retained. Set to your smallest experimental group size. Default 3.
#'
#' @return A fitted DESeqDataSet object.
build_dds <- function(counts, metadata, design,
                      min_count = 10, smallest_group_size = 3) {

  # Sanity checks
  if (!all(rownames(metadata) %in% colnames(counts))) {
    missing <- setdiff(rownames(metadata), colnames(counts))
    stop("Some metadata samples are not in the counts matrix: ",
         paste(missing, collapse = ", "))
  }

  # Subset and reorder counts to match metadata
  counts <- counts[, rownames(metadata)]
  stopifnot(all(colnames(counts) == rownames(metadata)))

  cat("Building DESeq2 object:\n")
  cat("  Samples:        ", ncol(counts), "\n")
  cat("  Initial genes:  ", nrow(counts), "\n")
  cat("  Design:         ", deparse(design), "\n")

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = design)

  # Pre-filter low-count genes
  keep <- rowSums(counts(dds) >= min_count) >= smallest_group_size
  dds <- dds[keep, ]
  cat("  Genes after filter (>=", min_count, "in >=", smallest_group_size,
      "samples):", nrow(dds), "\n")

  cat("Running DESeq() ...\n")
  dds <- DESeq(dds)

  cat("Available result names:\n")
  print(resultsNames(dds))

  dds
}


# --------------------------------------------------------------------------
# extract_contrast()
# --------------------------------------------------------------------------
#' Extract a clean results table from a fitted DESeq2 object
#'
#' Applies LFC shrinkage (apeglm if a coef name is given, ashr if a contrast
#' vector is given), joins gene annotations, computes per-group means, and
#' returns a tidy data frame sorted by adjusted p-value.
#'
#' @param dds A fitted DESeqDataSet (from build_dds()).
#' @param gene_info Gene annotation data.frame (from load_plasmidsaurus()$gene_info).
#' @param coef Character coefficient name from resultsNames(dds), e.g.
#'   "Exposure_Crystalline_vs_CTRL". Use this OR contrast, not both.
#' @param contrast Character vector of length 3: c(factor, numerator, denominator).
#'   Used when the contrast you want isn't a model coefficient. Triggers ashr shrinkage.
#' @param group_col Column in colData(dds) that defines groups for computing
#'   per-group normalized count means. Inferred from coef/contrast if NULL.
#' @param numerator Group label for the numerator of the contrast (for column naming).
#' @param denominator Group label for the denominator of the contrast.
#' @param alpha Significance level for the independent filtering step. Default 0.05.
#'
#' @return A data.frame with columns: gene_id, gene_name, gene_biotype,
#'   mean_<numerator>, mean_<denominator>, baseMean, log2FoldChange, lfcSE,
#'   pvalue, padj. Sorted by padj.
extract_contrast <- function(dds, gene_info,
                             coef = NULL, contrast = NULL,
                             group_col = NULL,
                             numerator = NULL, denominator = NULL,
                             alpha = 0.05) {

  # Decide shrinkage strategy
  if (!is.null(coef) && is.null(contrast)) {
    cat("Using coef-based shrinkage (apeglm) for:", coef, "\n")
    res_unshrunk <- results(dds, name = coef, alpha = alpha)
    res_shrunk   <- lfcShrink(dds, coef = coef, type = "apeglm")
    contrast_label <- coef
  } else if (!is.null(contrast) && is.null(coef)) {
    cat("Using contrast-based shrinkage (ashr) for:",
        paste(contrast, collapse = " "), "\n")
    res_unshrunk <- results(dds, contrast = contrast, alpha = alpha)
    res_shrunk   <- lfcShrink(dds, contrast = contrast, res = res_unshrunk,
                              type = "ashr")
    if (is.null(group_col))   group_col <- contrast[1]
    if (is.null(numerator))   numerator <- contrast[2]
    if (is.null(denominator)) denominator <- contrast[3]
    contrast_label <- paste(contrast, collapse = "_")
  } else {
    stop("Provide exactly one of `coef` or `contrast`.")
  }

  # Try to infer group column and labels from coef name (e.g. "Exposure_Crystalline_vs_CTRL")
  if (!is.null(coef)) {
    parts <- strsplit(coef, "_vs_")[[1]]
    if (length(parts) == 2) {
      if (is.null(denominator)) denominator <- parts[2]
      remaining <- strsplit(parts[1], "_")[[1]]
      if (is.null(group_col))   group_col <- remaining[1]
      if (is.null(numerator))   numerator <- paste(remaining[-1], collapse = "_")
    }
  }

  # Per-group normalized count means
  nc <- counts(dds, normalized = TRUE)
  if (!is.null(group_col) && group_col %in% colnames(colData(dds))) {
    grp <- colData(dds)[[group_col]]
    mean_num <- if (numerator %in% as.character(grp))
      rowMeans(nc[, grp == numerator, drop = FALSE]) else NA
    mean_den <- if (denominator %in% as.character(grp))
      rowMeans(nc[, grp == denominator, drop = FALSE]) else NA
  } else {
    mean_num <- NA; mean_den <- NA
  }

  # Build tidy table
  tbl <- as.data.frame(res_shrunk) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_info, by = "gene_id") %>%
    dplyr::mutate(
      !!paste0("mean_", numerator)   := mean_num[gene_id],
      !!paste0("mean_", denominator) := mean_den[gene_id]
    ) %>%
    dplyr::select(gene_id, gene_name, gene_biotype,
                  dplyr::any_of(c(paste0("mean_", numerator),
                                  paste0("mean_", denominator))),
                  baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
    dplyr::arrange(padj)

  attr(tbl, "contrast_label") <- contrast_label
  attr(tbl, "numerator")      <- numerator
  attr(tbl, "denominator")    <- denominator
  tbl
}


# --------------------------------------------------------------------------
# summarize_significance()
# --------------------------------------------------------------------------
#' Print summary counts at standard significance thresholds
#'
#' @param tbl Results table from extract_contrast().
#' @param padj_thresh Adjusted p-value threshold. Default 0.05.
#' @param lfc_thresh Absolute log2 fold change threshold. Default 1.
summarize_significance <- function(tbl, padj_thresh = 0.05, lfc_thresh = 1) {
  tbl_clean <- tbl %>% dplyr::filter(!is.na(padj))
  total <- nrow(tbl_clean)

  cat("=== Significance summary ===\n")
  cat("Genes tested (non-NA padj):", total, "\n\n")
  cat("padj <", padj_thresh, ":\n")
  cat("  total:    ", sum(tbl_clean$padj < padj_thresh), "\n")

  sig_strict <- tbl_clean %>%
    dplyr::filter(padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh)
  sig_loose <- tbl_clean %>%
    dplyr::filter(padj < padj_thresh, abs(log2FoldChange) >= 0.585)

  cat("\npadj <", padj_thresh, "AND |LFC| >=", lfc_thresh,
      "(", round(2^lfc_thresh, 2), "-fold):\n")
  cat("  total:    ", nrow(sig_strict), "\n")
  cat("  upregulated:  ", sum(sig_strict$log2FoldChange > 0), "\n")
  cat("  downregulated:", sum(sig_strict$log2FoldChange < 0), "\n")

  cat("\npadj <", padj_thresh, "AND |LFC| >= 0.585 (1.5-fold, more lenient):\n")
  cat("  total:    ", nrow(sig_loose), "\n")

  invisible(list(strict = sig_strict, loose = sig_loose))
}


# --------------------------------------------------------------------------
# plot_pca_custom()
# --------------------------------------------------------------------------
#' PCA plot colored by any metadata column(s)
#'
#' @param rld An rlog- or vst-transformed DESeqTransform object.
#' @param color_by Metadata column to color points by.
#' @param shape_by Optional metadata column for point shape.
#' @param ntop Number of most variable genes to use. Default 5000.
#' @param label Logical, whether to label points with sample names. Default TRUE.
#' @param title Plot title.
#'
#' @return A ggplot object.
plot_pca_custom <- function(rld, color_by, shape_by = NULL,
                            ntop = 5000, label = TRUE,
                            title = "PCA on rlog-transformed counts") {
  intgroup <- c(color_by, shape_by)
  pca_data <- plotPCA(rld, intgroup = intgroup, ntop = ntop, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

  aes_args <- list(x = quote(PC1), y = quote(PC2),
                   color = as.name(color_by), label = quote(name))
  if (!is.null(shape_by)) aes_args$shape <- as.name(shape_by)

  p <- ggplot(pca_data, do.call(aes, aes_args)) +
    geom_point(size = 4) +
    xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    theme_bw(base_size = 14) +
    ggtitle(title)

  if (label) {
    p <- p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE)
  }
  p
}


# --------------------------------------------------------------------------
# plot_sample_distance()
# --------------------------------------------------------------------------
#' Sample-to-sample Euclidean distance heatmap
#'
#' @param rld An rlog- or vst-transformed DESeqTransform object.
#' @param label_cols Character vector of metadata columns to use in row labels.
#'
#' @return Invisibly, the pheatmap object.
plot_sample_distance <- function(rld, label_cols = NULL) {
  d <- dist(t(assay(rld)))
  m <- as.matrix(d)

  if (!is.null(label_cols) && all(label_cols %in% colnames(colData(rld)))) {
    labs <- apply(as.data.frame(colData(rld)[, label_cols, drop = FALSE]),
                  1, paste, collapse = " | ")
    rownames(m) <- paste(labs, rownames(m), sep = " - ")
  }
  colnames(m) <- NULL
  cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

  pheatmap(m,
           clustering_distance_rows = d,
           clustering_distance_cols = d,
           col = cols,
           main = "Sample-to-sample distance (rlog)")
}


# --------------------------------------------------------------------------
# plot_volcano_custom()
# --------------------------------------------------------------------------
#' Volcano plot with sensible defaults using EnhancedVolcano
#'
#' @param tbl Results table from extract_contrast().
#' @param padj_thresh Significance threshold. Default 0.05.
#' @param lfc_thresh Log2 fold change threshold. Default 1.
#' @param title Plot title (uses contrast_label attribute if NULL).
#' @param top_n_label Top N genes to label by padj. Default 20.
#'
#' @return An EnhancedVolcano ggplot.
plot_volcano_custom <- function(tbl, padj_thresh = 0.05, lfc_thresh = 1,
                                title = NULL, top_n_label = 20) {
  if (is.null(title)) {
    title <- attr(tbl, "contrast_label")
    if (is.null(title)) title <- "Volcano plot"
  }

  # Pick top N labeled genes by padj
  tbl_with_labels <- tbl
  tbl_with_labels$display_label <- ""
  top_idx <- order(tbl$padj, na.last = TRUE)[seq_len(min(top_n_label, nrow(tbl)))]
  tbl_with_labels$display_label[top_idx] <- tbl$gene_name[top_idx]

  EnhancedVolcano::EnhancedVolcano(
    tbl_with_labels,
    lab = tbl_with_labels$display_label,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = padj_thresh,
    FCcutoff = lfc_thresh,
    title = title,
    subtitle = paste0("Cutoffs: padj < ", padj_thresh,
                      ", |LFC| >= ", lfc_thresh),
    pointSize = 2,
    labSize = 4,
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    colAlpha = 0.7,
    legendPosition = "right"
  )
}


# --------------------------------------------------------------------------
# write_results_excel()
# --------------------------------------------------------------------------
#' Write a results table to a formatted Excel file
#'
#' @param tbl Results table from extract_contrast().
#' @param path Output file path (should end in .xlsx).
#' @param sheet_name Sheet name. Default "DEGs".
write_results_excel <- function(tbl, path, sheet_name = "DEGs") {
  openxlsx::write.xlsx(tbl, file = path,
                       sheetName = sheet_name,
                       asTable = TRUE,
                       tableStyle = "TableStyleMedium2",
                       firstRow = TRUE,
                       headerStyle = openxlsx::createStyle(
                         textDecoration = "bold",
                         fgFill = "#4F81BD",
                         fontColour = "white"))
  cat("Wrote", path, "(", nrow(tbl), "genes)\n")
}


# --------------------------------------------------------------------------
# write_outputs_for_enrichment()
# --------------------------------------------------------------------------
#' Save standard outputs needed by functional enrichment tools
#'
#' Writes:
#'   - <prefix>_full_results.csv: full results table (also as .xlsx)
#'   - <prefix>_significant.csv:  significant genes only
#'   - <prefix>_ranked.csv:       all genes ranked by -log10(p) * sign(LFC) for GSEA
#'
#' @param tbl Results table from extract_contrast().
#' @param outdir Output directory.
#' @param prefix File prefix.
#' @param padj_thresh Significance threshold. Default 0.05.
#' @param lfc_thresh LFC threshold for "significant" list. Default 1.
write_outputs_for_enrichment <- function(tbl, outdir, prefix,
                                          padj_thresh = 0.05, lfc_thresh = 1) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Full results
  write.csv(tbl, file.path(outdir, paste0(prefix, "_full_results.csv")),
            row.names = FALSE)
  write_results_excel(tbl, file.path(outdir, paste0(prefix, "_full_results.xlsx")))

  # Significant only
  sig <- tbl %>%
    dplyr::filter(!is.na(padj), padj < padj_thresh,
                  abs(log2FoldChange) >= lfc_thresh)
  write.csv(sig, file.path(outdir, paste0(prefix, "_significant.csv")),
            row.names = FALSE)

  # Ranked for GSEA
  ranked <- tbl %>%
    dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
    dplyr::mutate(rank_metric = -log10(pvalue) * sign(log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(rank_metric)) %>%
    dplyr::select(gene_id, gene_name, log2FoldChange, padj, rank_metric)
  write.csv(ranked, file.path(outdir, paste0(prefix, "_ranked.csv")),
            row.names = FALSE)

  cat("Saved enrichment-ready outputs to", outdir, "with prefix", prefix, "\n")
  invisible(list(full = tbl, significant = sig, ranked = ranked))
}


# ============================================================================
# QC FUNCTIONS
# ============================================================================
# Run BEFORE the differential expression analysis. The goal of this section is
# to surface anything weird about your samples (low depth, unusual biotype
# composition, mislabeling, low correlation with replicates) so you can decide
# whether to drop samples before they distort downstream results.
#
# Typical use from a notebook:
#   qc <- run_qc_report(
#     counts        = ps$counts,
#     gene_info     = ps$gene_info,
#     metadata      = metadata,
#     mapping_stats = "raw_data/SS8VZG-mapping-stats-reads.csv",
#     output_dir    = "qc"
#   )
#   # Then look at qc$flagged_samples and decide what to do
# ============================================================================


# --------------------------------------------------------------------------
# load_mapping_stats()
# --------------------------------------------------------------------------
#' Load Plasmidsaurus mapping stats CSV
#'
#' @param path Path to the SS8VZG-mapping-stats-reads.csv file (or similar).
#' @return data.frame with columns: sample_label, Uniquely Mapped, Multi-mapped,
#'   Unmapped, Total, pct_unique, pct_multi, pct_unmapped
load_mapping_stats <- function(path) {
  if (!file.exists(path)) {
    warning("Mapping stats file not found: ", path,
            " - mapping QC will be skipped.")
    return(NULL)
  }
  ms <- read.csv(path, check.names = FALSE)
  colnames(ms)[1] <- "sample_label"
  ms$Total          <- ms$`Uniquely Mapped` + ms$`Multi-mapped` + ms$Unmapped
  ms$pct_unique     <- 100 * ms$`Uniquely Mapped` / ms$Total
  ms$pct_multi      <- 100 * ms$`Multi-mapped`   / ms$Total
  ms$pct_unmapped   <- 100 * ms$Unmapped         / ms$Total
  ms
}


# --------------------------------------------------------------------------
# compute_qc_metrics()
# --------------------------------------------------------------------------
#' Compute per-sample QC metrics from counts and gene annotations
#'
#' @param counts Integer counts matrix (genes x samples), rownames = gene_id.
#' @param gene_info data.frame with columns gene_id, gene_name, gene_biotype.
#' @param metadata data.frame with row per sample. Optional but useful for grouping.
#' @return data.frame with one row per sample. Columns include:
#'   sample_id, library_size, n_genes_detected_1, n_genes_detected_5,
#'   n_genes_detected_10, pct_protein_coding, pct_mt, pct_ribosomal,
#'   top_gene_pct (% of reads in the single most-expressed gene)
compute_qc_metrics <- function(counts, gene_info, metadata = NULL) {

  # Identify gene categories using gene_info and gene_name patterns
  is_protein <- gene_info$gene_biotype == "protein_coding"
  # Mitochondrial genes in mouse Ensembl have prefix "mt-" in gene_name
  is_mito    <- grepl("^mt-", gene_info$gene_name, ignore.case = TRUE)
  # Ribosomal protein genes (Rps*, Rpl*, Mrps*, Mrpl*)
  is_ribo    <- grepl("^Rp[sl][0-9]|^Mrp[sl][0-9]",
                      gene_info$gene_name, ignore.case = TRUE)

  # Match gene_info ordering to counts ordering
  m <- match(rownames(counts), gene_info$gene_id)
  is_protein <- is_protein[m]
  is_mito    <- is_mito[m]
  is_ribo    <- is_ribo[m]

  out <- data.frame(
    sample_id            = colnames(counts),
    library_size         = colSums(counts),
    n_genes_detected_1   = colSums(counts >= 1),
    n_genes_detected_5   = colSums(counts >= 5),
    n_genes_detected_10  = colSums(counts >= 10),
    pct_protein_coding   = 100 * colSums(counts[is_protein, , drop = FALSE]) /
                                  pmax(colSums(counts), 1),
    pct_mt               = 100 * colSums(counts[is_mito, , drop = FALSE], na.rm = TRUE) /
                                  pmax(colSums(counts), 1),
    pct_ribosomal        = 100 * colSums(counts[is_ribo, , drop = FALSE], na.rm = TRUE) /
                                  pmax(colSums(counts), 1),
    top_gene_pct         = 100 * apply(counts, 2, max) / pmax(colSums(counts), 1),
    stringsAsFactors = FALSE
  )

  # Attach metadata columns (Sex, Exposure, Tissue, etc.) if available
  if (!is.null(metadata)) {
    md <- metadata[match(out$sample_id, rownames(metadata)), , drop = FALSE]
    # Drop sample_id from metadata to avoid duplicate column after cbind
    md$sample_id <- NULL
    out <- cbind(out, md)
  }

  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------
# check_sex_markers()
# --------------------------------------------------------------------------
#' Check expression of sex-linked genes to verify sex assignment
#'
#' Females express XIST/Xist (X-inactive specific transcript). Males express
#' Y-linked genes. Mismatches between expression and the declared sex in
#' metadata indicate a labeling error.
#'
#' @param counts Counts matrix.
#' @param gene_info data.frame with gene_id, gene_name.
#' @param metadata data.frame with a Sex column (values "F"/"M" or similar).
#' @param species Species name. "mouse" / "Mus musculus" or "human" / "Homo sapiens".
#'   Default "mouse" for backward compatibility. Determines whether to use
#'   lowercase mouse gene symbols (Xist, Ddx3y) or uppercase human (XIST, DDX3Y).
#' @param normalize Logical, whether to library-size-normalize the counts before
#'   comparison. Default TRUE.
#' @return data.frame: per-sample expression of each sex marker, plus declared sex.
check_sex_markers <- function(counts, gene_info, metadata,
                              species = "mouse", normalize = TRUE) {

  species_l <- tolower(species)
  if (species_l %in% c("mouse", "mus musculus", "mm")) {
    markers_F <- c("Xist", "Tsix")
    markers_M <- c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d")
    primary_F <- "Xist"
  } else if (species_l %in% c("human", "homo sapiens", "hs")) {
    markers_F <- c("XIST")
    markers_M <- c("RPS4Y1", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y", "NLGN4Y")
    primary_F <- "XIST"
  } else {
    stop("Unknown species: ", species,
         ". Use 'mouse'/'Mus musculus' or 'human'/'Homo sapiens'.")
  }
  all_markers <- c(markers_F, markers_M)

  # Find their gene_ids
  ids <- gene_info$gene_id[match(all_markers, gene_info$gene_name)]
  found <- !is.na(ids)
  if (!any(found)) {
    warning("No sex markers found in gene_info. Check organism/annotation.")
    return(NULL)
  }
  marker_counts <- counts[ids[found], , drop = FALSE]
  rownames(marker_counts) <- all_markers[found]

  if (normalize) {
    lib <- colSums(counts)
    marker_counts <- sweep(marker_counts, 2, lib / 1e6, "/")  # CPM
  }

  out <- as.data.frame(t(marker_counts))
  out$sample_id <- colnames(marker_counts)
  out$declared_sex <- metadata$Sex[match(out$sample_id, rownames(metadata))]
  out <- out[, c("sample_id", "declared_sex", rownames(marker_counts))]

  # Heuristic flag: F samples should have high X-linked and low Y-genes; vice versa for M
  female_score <- if (primary_F %in% rownames(marker_counts))
    out[[primary_F]] else 0
  m_present <- intersect(markers_M, colnames(out))
  male_score <- if (length(m_present) > 0)
    rowMeans(out[, m_present, drop = FALSE], na.rm = TRUE) else rep(0, nrow(out))

  out$predicted_sex <- ifelse(female_score > male_score, "F", "M")
  out$sex_consistent <- out$declared_sex == out$predicted_sex

  rownames(out) <- NULL
  out
}


# --------------------------------------------------------------------------
# plot_qc_summary()
# --------------------------------------------------------------------------
#' Multi-panel QC figure
#'
#' @param qc_df Output of compute_qc_metrics().
#' @param mapping_stats Output of load_mapping_stats(), or NULL.
#' @param color_by Metadata column to color bars by. Default "Exposure" if present.
#' @return A patchwork ggplot object (or just one ggplot if patchwork unavailable).
plot_qc_summary <- function(qc_df, mapping_stats = NULL, color_by = NULL) {

  if (is.null(color_by)) {
    color_by <- if ("Exposure" %in% colnames(qc_df)) "Exposure" else NULL
  }

  base_aes <- function(y_col) {
    if (!is.null(color_by) && color_by %in% colnames(qc_df)) {
      aes(x = reorder(sample_id, library_size), y = .data[[y_col]],
          fill = .data[[color_by]])
    } else {
      aes(x = reorder(sample_id, library_size), y = .data[[y_col]])
    }
  }

  theme_qc <- theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          plot.title  = element_text(size = 11),
          legend.position = "bottom")

  p_lib <- ggplot(qc_df, base_aes("library_size")) +
    geom_col() +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Library size (total reads)", x = NULL, y = "Reads") +
    theme_qc

  p_det <- ggplot(qc_df, base_aes("n_genes_detected_10")) +
    geom_col() +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Genes detected (\u2265 10 reads)", x = NULL, y = "Genes") +
    theme_qc

  p_pc <- ggplot(qc_df, base_aes("pct_protein_coding")) +
    geom_col() +
    geom_hline(yintercept = 70, linetype = "dashed", color = "grey40") +
    labs(title = "% protein-coding (>=70% expected)",
         x = NULL, y = "% reads in protein-coding") +
    theme_qc

  p_mt <- ggplot(qc_df, base_aes("pct_mt")) +
    geom_col() +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    labs(title = "% mitochondrial (>10% suggests degradation)",
         x = NULL, y = "% reads in mt-genes") +
    theme_qc

  p_top <- ggplot(qc_df, base_aes("top_gene_pct")) +
    geom_col() +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    labs(title = "% reads in single top gene (>10% = dominant gene)",
         x = NULL, y = "%") +
    theme_qc

  plots <- list(p_lib, p_det, p_pc, p_mt, p_top)

  if (!is.null(mapping_stats)) {
    # Match mapping stats to qc_df sample_ids
    # The mapping stats use the original sample labels; qc_df uses stripped IDs.
    # Plasmidsaurus's CSV may inconsistently use either '-' or '_' as the separator
    # (and even mix them within one file), so we extract the leading numeric ID
    # by stopping at the first separator of either type.
    ms <- mapping_stats
    ms$sample_id <- as.character(sub("[-_].*$", "", ms$sample_label))

    if (any(ms$sample_id %in% qc_df$sample_id)) {
      ms_match <- ms[match(qc_df$sample_id, ms$sample_id), ]
      qc_df$pct_unique   <- ms_match$pct_unique
      qc_df$pct_unmapped <- ms_match$pct_unmapped

      p_uniq <- ggplot(qc_df, base_aes("pct_unique")) +
        geom_col() +
        geom_hline(yintercept = 70, linetype = "dashed", color = "grey40") +
        labs(title = "% uniquely mapped (>=70% expected)",
             x = NULL, y = "%") +
        theme_qc

      p_unmap <- ggplot(qc_df, base_aes("pct_unmapped")) +
        geom_col() +
        geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
        labs(title = "% unmapped (>10% suggests problems)",
             x = NULL, y = "%") +
        theme_qc

      plots <- c(plots, list(p_uniq, p_unmap))
    }
  }

  # Combine with patchwork if available, else return list
  if (requireNamespace("patchwork", quietly = TRUE)) {
    Reduce(`/`, plots[-1], plots[[1]])
  } else {
    cat("Install 'patchwork' for combined plots. Returning list.\n")
    plots
  }
}


# --------------------------------------------------------------------------
# plot_sample_correlation()
# --------------------------------------------------------------------------
#' Pearson correlation heatmap between samples (on log2-transformed counts)
#'
#' Replicates within a group should correlate >= 0.9. Lower than ~0.85 between
#' biological replicates is a warning sign.
#'
#' @param counts Integer counts matrix.
#' @param metadata Optional data.frame for annotation. Useful column names are
#'   auto-detected (Group/Exposure, Sex, Tissue/tissue). Use `ann_cols` to
#'   override which columns to show.
#' @param ann_cols Character vector of metadata column names to show as annotations.
#'   If NULL, auto-detects sensible columns.
#' @param min_count Genes must have >= this in at least one sample to be included. Default 10.
#' @return Invisibly, a list with cor_matrix (the correlation matrix) and the pheatmap object.
plot_sample_correlation <- function(counts, metadata = NULL,
                                    ann_cols = NULL, min_count = 10) {
  keep <- rowSums(counts >= min_count) > 0
  filtered <- counts[keep, , drop = FALSE]
  log_counts <- log2(filtered + 1)
  cor_mat <- cor(log_counts, method = "pearson")

  ann <- NULL
  if (!is.null(metadata)) {
    if (is.null(ann_cols)) {
      # Auto-detect sensible annotation columns
      candidate <- c("Group", "Exposure", "Sex", "Tissue", "tissue",
                     "Condition", "Genotype", "genotype")
      ann_cols <- intersect(candidate, colnames(metadata))
    } else {
      ann_cols <- intersect(ann_cols, colnames(metadata))
    }

    if (length(ann_cols) > 0) {
      ann_raw <- metadata[colnames(cor_mat), ann_cols, drop = FALSE]
      # Drop columns that are all NA or all empty (pheatmap will choke on these)
      keep_ann <- sapply(ann_raw, function(x) {
        x_chr <- as.character(x)
        !all(is.na(x_chr) | x_chr == "")
      })
      ann_raw <- ann_raw[, keep_ann, drop = FALSE]
      if (ncol(ann_raw) > 0) {
        ann <- ann_raw
      }
    }
  }

  ph <- pheatmap(cor_mat,
                 annotation_col = ann,
                 annotation_row = ann,
                 main = "Sample-sample Pearson correlation (log2 counts)",
                 display_numbers = TRUE, number_format = "%.2f",
                 fontsize_number = 6,
                 color = colorRampPalette(c("white", "#fee0d2", "#fc9272",
                                            "#de2d26"))(100))
  invisible(list(cor_matrix = cor_mat, plot = ph))
}


# --------------------------------------------------------------------------
# flag_qc_outliers()
# --------------------------------------------------------------------------
#' Apply thresholds to qc metrics and flag samples for review
#'
#' Default thresholds reflect generally accepted RNA-seq QC norms; tune as
#' needed for your data. Each threshold can be overridden.
#'
#' @param qc_df Output of compute_qc_metrics().
#' @param min_library Minimum library size. Default 5e6 (5M reads).
#' @param min_genes Minimum genes detected at >= 10 reads. Default 12000.
#' @param min_pct_protein_coding Minimum % protein-coding. Default 60.
#' @param max_pct_mt Maximum % mitochondrial. Default 15.
#' @param max_top_gene_pct Maximum % reads in the single top gene. Default 15.
#' @return data.frame with all qc_df columns plus per-threshold pass flags
#'   and a final pass_qc column. Samples with pass_qc = FALSE failed >=1 check.
flag_qc_outliers <- function(qc_df,
                              min_library            = 5e6,
                              min_genes              = 12000,
                              min_pct_protein_coding = 60,
                              max_pct_mt             = 15,
                              max_top_gene_pct       = 15) {

  out <- qc_df
  out$pass_library      <- out$library_size >= min_library
  out$pass_detection    <- out$n_genes_detected_10 >= min_genes
  out$pass_protein      <- out$pct_protein_coding >= min_pct_protein_coding
  out$pass_mt           <- out$pct_mt <= max_pct_mt
  out$pass_top_gene     <- out$top_gene_pct <= max_top_gene_pct

  pass_cols <- c("pass_library", "pass_detection", "pass_protein",
                 "pass_mt", "pass_top_gene")
  out$pass_qc <- apply(out[, pass_cols], 1, all)

  cat("Thresholds applied:\n")
  cat("  library_size         >=", format(min_library, big.mark = ","), "\n")
  cat("  n_genes_detected_10  >=", min_genes, "\n")
  cat("  pct_protein_coding   >=", min_pct_protein_coding, "%\n")
  cat("  pct_mt               <=", max_pct_mt, "%\n")
  cat("  top_gene_pct         <=", max_top_gene_pct, "%\n\n")

  cat("Samples flagged for review (failing one or more checks):\n")
  flagged <- out[!out$pass_qc, ]
  if (nrow(flagged) == 0) {
    cat("  (none)\n")
  } else {
    for (i in seq_len(nrow(flagged))) {
      reasons <- c()
      if (!flagged$pass_library[i])   reasons <- c(reasons, sprintf("low library (%s)",
                                                    format(flagged$library_size[i], big.mark = ",")))
      if (!flagged$pass_detection[i]) reasons <- c(reasons, sprintf("low detection (%d)",
                                                    flagged$n_genes_detected_10[i]))
      if (!flagged$pass_protein[i])   reasons <- c(reasons, sprintf("low %% protein (%.1f%%)",
                                                    flagged$pct_protein_coding[i]))
      if (!flagged$pass_mt[i])        reasons <- c(reasons, sprintf("high %% mt (%.1f%%)",
                                                    flagged$pct_mt[i]))
      if (!flagged$pass_top_gene[i])  reasons <- c(reasons, sprintf("dominant gene (%.1f%%)",
                                                    flagged$top_gene_pct[i]))
      cat("  ", flagged$sample_id[i], ":", paste(reasons, collapse = "; "), "\n")
    }
  }
  cat("\nNote: flagging is advisory. Review the QC plots before deciding.\n")

  out
}


# --------------------------------------------------------------------------
# run_qc_report()
# --------------------------------------------------------------------------
#' One-call wrapper that runs the full QC report and saves outputs
#'
#' @param counts Counts matrix from load_plasmidsaurus()$counts.
#' @param gene_info Gene info from load_plasmidsaurus()$gene_info.
#' @param metadata data.frame with sample metadata. rownames = sample_id.
#' @param mapping_stats Path to mapping stats CSV, or NULL.
#' @param output_dir Directory to write QC outputs.
#' @param thresholds Named list to override default flag_qc_outliers() thresholds.
#' @param species Species for sex marker check ("mouse" or "human"). Default "mouse".
#' @return list with elements: metrics, sex_check, flagged_samples
run_qc_report <- function(counts, gene_info, metadata,
                          mapping_stats = NULL,
                          output_dir = "qc",
                          thresholds = list(),
                          species = "mouse") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. Compute basic metrics
  cat("=== Computing QC metrics ===\n")
  metrics <- compute_qc_metrics(counts, gene_info, metadata)

  # 2. Mapping stats (optional, Plasmidsaurus-specific)
  ms <- NULL
  if (!is.null(mapping_stats)) {
    cat("=== Loading mapping stats ===\n")
    ms <- load_mapping_stats(mapping_stats)
  }

  # 3. Sex check (skip if there's no Sex column or all values are NA)
  sex_check <- NULL
  has_sex <- "Sex" %in% colnames(metadata) &&
             !all(is.na(metadata$Sex)) &&
             !all(metadata$Sex == "")
  if (has_sex) {
    cat("=== Checking sex markers ===\n")
    sex_check <- check_sex_markers(counts, gene_info, metadata, species = species)
    if (!is.null(sex_check)) {
      mismatches <- sex_check[!sex_check$sex_consistent &
                              !is.na(sex_check$sex_consistent), ]
      if (nrow(mismatches) > 0) {
        cat("WARNING: Sex marker mismatches detected!\n")
        print(mismatches[, c("sample_id", "declared_sex", "predicted_sex")])
      } else {
        cat("All", nrow(sex_check), "samples consistent with declared sex.\n")
      }
    }
  } else {
    cat("=== Sex check skipped (no Sex column or all NA) ===\n")
  }

  # 4. Flag outliers
  cat("\n=== Flagging samples ===\n")
  flag_args <- c(list(qc_df = metrics), thresholds)
  flagged <- do.call(flag_qc_outliers, flag_args)

  # 5. Save tables
  write.csv(flagged, file.path(output_dir, "qc_metrics.csv"), row.names = FALSE)
  if (!is.null(sex_check)) {
    write.csv(sex_check, file.path(output_dir, "qc_sex_check.csv"), row.names = FALSE)
  }
  if (!is.null(ms)) {
    write.csv(ms, file.path(output_dir, "qc_mapping_stats.csv"), row.names = FALSE)
  }

  cat("\nQC outputs saved to:", output_dir, "\n")
  cat("Files:\n")
  cat("  qc_metrics.csv     - Per-sample metrics with pass/fail flags\n")
  if (!is.null(sex_check)) cat("  qc_sex_check.csv   - Sex marker expression and consistency\n")
  if (!is.null(ms))        cat("  qc_mapping_stats.csv - Mapping statistics\n")

  invisible(list(metrics = flagged,
                 sex_check = sex_check,
                 mapping_stats = ms,
                 flagged_samples = flagged$sample_id[!flagged$pass_qc]))
}
