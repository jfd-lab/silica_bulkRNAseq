# ============================================================================
# load_featurecounts.R
# ============================================================================
#
# Drop-in replacement for load_plasmidsaurus() that reads the merged
# featureCounts matrix produced by your Bouchet pipeline (5_merge_counts.r ->
# rawcounts_GRCm39.txt) instead of the Plasmidsaurus expression-matrix TSV.
#
# It returns the SAME structure load_plasmidsaurus() returns, so every
# downstream function in pipeline_utils.R (compute_qc_metrics, check_sex_markers,
# build_dds, extract_contrast, etc.) works unchanged:
#
#   list(counts = <integer matrix, rownames = gene_id>,
#        gene_info = data.frame(gene_id, gene_name, gene_biotype),
#        sample_ids = <chr>)
#
# KEY DIFFERENCE from Plasmidsaurus:
#   featureCounts here counted on gene_name (GTF.attrType="gene_name"), so the
#   GeneID column is a gene SYMBOL, not an Ensembl ID. We therefore set BOTH
#   gene_id and gene_name to that symbol. Downstream code keys on gene_id for
#   joins/rownames and on gene_name for mito/sex/ribo pattern matching -- both
#   work because they're the same symbol string here.
#
#   featureCounts output has NO biotype column. compute_qc_metrics() uses
#   gene_biotype for pct_protein_coding. Two modes:
#     - gtf_path = NULL (default): gene_biotype = NA. pct_protein_coding will be
#       NA/0; either ignore that one QC metric or relax its threshold (see notes
#       at bottom).
#     - gtf_path = "<your GENCODE gtf>": we parse gene_type from the GTF and fill
#       gene_biotype properly, so ALL QC metrics work exactly as before.
#
# Usage in your notebook (replaces the load_plasmidsaurus call):
#   source("pipeline_utils.R")
#   source("load_featurecounts.R")
#   ps <- load_featurecounts(
#     counts_path = "/nfs/roberts/project/pi_epm38/jfd42/silica/rawcounts_GRCm39.txt",
#     gtf_path    = "/nfs/roberts/project/pi_epm38/jfd42/GENCODE/GRCm39/gencode.vM39.primary_assembly.annotation.gtf"
#   )
#   # ps$counts, ps$gene_info, ps$sample_ids  -- identical shape to before
# ============================================================================

load_featurecounts <- function(counts_path,
                               gtf_path = NULL,
                               strip_bam_suffix = TRUE,
                               strip_prefix = NULL) {

  if (!file.exists(counts_path)) {
    stop("File not found: ", counts_path)
  }

  cat("Reading", counts_path, "...\n")
  raw <- read.delim(counts_path, sep = "\t", check.names = FALSE,
                    stringsAsFactors = FALSE)

  # featureCounts annotation columns (from Rsubread featureCounts$annotation)
  annot_cols <- c("GeneID", "Chr", "Start", "End", "Strand", "Length")
  missing <- setdiff(annot_cols, colnames(raw))
  if (length(missing) > 0) {
    stop("Expected featureCounts columns missing: ", paste(missing, collapse = ", "),
         "\n(Is this the merged rawcounts_*.txt from 5_merge_counts.r?)")
  }

  # Count columns = everything that is not an annotation column
  count_cols <- setdiff(colnames(raw), annot_cols)
  if (length(count_cols) == 0) {
    stop("No count columns found (only annotation columns present).")
  }

  # Build counts matrix, rownames = GeneID (the gene symbol)
  counts <- as.matrix(raw[, count_cols, drop = FALSE])
  rownames(counts) <- raw$GeneID

  # Clean sample names ---------------------------------------------------------
  # featureCounts headers can carry the BAM filename, e.g.
  #   "F1C2Aligned.sortedByCoord.out.bam" -> "F1C2"
  if (strip_bam_suffix) {
    colnames(counts) <- sub("Aligned.*$", "", colnames(counts))
    colnames(counts) <- sub("\\.bam$", "", colnames(counts))
  }
  if (!is.null(strip_prefix)) {
    colnames(counts) <- sub(strip_prefix, "", colnames(counts))
  }

  # featureCounts is already integer, but enforce it for DESeq2
  storage.mode(counts) <- "integer"

  # Build gene_info to match the load_plasmidsaurus() contract -----------------
  gene_info <- data.frame(
    gene_id      = raw$GeneID,   # symbol (that's what we counted on)
    gene_name    = raw$GeneID,   # same symbol; downstream mito/sex/ribo regexes use this
    gene_biotype = NA_character_,
    stringsAsFactors = FALSE
  )

  # Optionally enrich biotype from the GTF (gene_type attribute) ---------------
  if (!is.null(gtf_path)) {
    if (!file.exists(gtf_path)) {
      warning("gtf_path not found: ", gtf_path,
              " -- gene_biotype will remain NA.")
    } else {
      cat("Parsing gene_type from GTF for biotype annotation...\n")
      bt <- .parse_gene_biotype_from_gtf(gtf_path)
      # Map by gene_name (symbol) since that's our key
      idx <- match(gene_info$gene_name, bt$gene_name)
      gene_info$gene_biotype <- bt$gene_biotype[idx]
      n_annotated <- sum(!is.na(gene_info$gene_biotype))
      cat("  Biotype assigned for", n_annotated, "of", nrow(gene_info), "genes\n")
    }
  } else {
    cat("NOTE: no gtf_path given -> gene_biotype = NA.\n")
    cat("      pct_protein_coding QC metric will be NA. See notes in this file.\n")
  }

  # Report (mirrors load_plasmidsaurus output) --------------------------------
  cat("  Genes:  ", nrow(counts), "\n")
  cat("  Samples:", ncol(counts), "\n")
  cat("  Sample names:", paste(colnames(counts), collapse = ", "), "\n")
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
# Internal: parse gene_name -> gene_type from a GENCODE GTF
# --------------------------------------------------------------------------
# Reads only the 'gene' feature lines and extracts gene_name + gene_type.
# Lightweight (no GTF-parsing packages needed); one pass with regex.
.parse_gene_biotype_from_gtf <- function(gtf_path) {
  con <- file(gtf_path, "r")
  on.exit(close(con))

  names_v  <- character(0)
  types_v  <- character(0)
  chunk_size <- 100000

  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) break
    # Keep only 'gene' feature rows (3rd tab field == "gene")
    lines <- lines[!startsWith(lines, "#")]
    if (length(lines) == 0) next
    # Fast field check: split is expensive; filter with grepl on "\tgene\t"
    gene_lines <- lines[grepl("\tgene\t", lines, fixed = TRUE)]
    if (length(gene_lines) == 0) next

    gn <- sub('.*gene_name "([^"]+)".*', "\\1", gene_lines)
    gt <- sub('.*gene_type "([^"]+)".*', "\\1", gene_lines)
    # Rows where extraction failed will equal the whole line; drop them
    ok <- gn != gene_lines & gt != gene_lines
    names_v <- c(names_v, gn[ok])
    types_v <- c(types_v, gt[ok])
  }

  df <- data.frame(gene_name = names_v, gene_biotype = types_v,
                   stringsAsFactors = FALSE)
  # De-duplicate (a symbol should map to one gene_type; keep first)
  df[!duplicated(df$gene_name), , drop = FALSE]
}


# ============================================================================
# NOTES / GOTCHAS
# ============================================================================
#
# 1. BIOTYPE: strongly recommend passing gtf_path so pct_protein_coding works.
#    Without it, in run_qc_report()/flag_qc_outliers() the pass_protein check
#    will fail for all samples (NA >= 60 is NA/FALSE). If you deliberately skip
#    the GTF, disable that one check by passing a permissive threshold, e.g.:
#       run_qc_report(..., thresholds = list(min_pct_protein_coding = 0))
#
# 2. GENE ID vs SYMBOL: because counting was on gene_name, a handful of Ensembl
#    genes that share a symbol (or have empty symbols) may have been collapsed
#    by featureCounts. This is inherent to Xiting's gene_name counting choice
#    (kept for methodological consistency). Your gene count (~78k) reflects that.
#    Duplicated symbols, if any, become duplicated rownames in `counts`; DESeq2
#    tolerates this but if you want strictly unique rownames, tell me and I'll
#    add a de-dup/aggregation step.
#
# 3. SEX MARKERS + MITO: these work unchanged. Mouse GENCODE mito genes are
#    "mt-Nd1" etc. (match ^mt-), and Xist/Ddx3y/etc. are present as symbols.
#
# 4. METADATA: unchanged. Your metadata.csv rownames must match the cleaned
#    sample names printed above (F1C2, CTRL_M1, ...). If they don't line up,
#    that's a name-cleaning issue, not a loader issue -- check the printout.
# ============================================================================
