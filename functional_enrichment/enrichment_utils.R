# ============================================================================
# enrichment_utils.R
# ============================================================================
#
# Reusable helper functions for functional enrichment analysis (ORA + GSEA)
# using clusterProfiler and MSigDB gene sets.
#
# Source this file from any analysis notebook with:
#   source("~/project_pi_epm38/jfd42/pipelines/functional_enrichment/enrichment_utils.R")
#
# Functions provided (in order of typical use):
#
#   --- Gene set loading ---
#   load_geneset_collection()  - Load one MSigDB collection (GO:BP, Hallmark, etc.)
#   load_all_collections()     - Load multiple collections at once
#
#   --- Input preparation ---
#   prepare_ora_input()        - From DESeq2 results -> ORA-ready gene lists
#   prepare_gsea_input()       - From DESeq2 results -> ranked named vector
#
#   --- Running enrichment ---
#   run_ora()                  - Wrapper around enricher() with our defaults
#   run_ora_all_collections()  - Run ORA across many collections + directions in one call
#   run_gsea()                 - Wrapper around GSEA() with our defaults
#   run_gsea_all_collections() - Run GSEA across many collections in one call
#
#   --- Plotting ---
#   plot_ora_dotplot()         - Standard dotplot
#   plot_ora_emap()            - Enrichment map (clusters related terms)
#   plot_gsea_dotplot()        - GSEA dotplot
#   plot_gsea_top_n()          - Grid of top N GSEA enrichment plots
#   compare_ora_dotplot()      - Side-by-side dotplot comparing multiple contrasts
#   plot_gsea_horizontal_bars()- Publication-style bidirectional bar plot (de Man style)
#
#   --- Output ---
#   save_enrichment_results()  - Save ORA/GSEA results to CSV + Excel with metadata
#   save_all_ora_results()     - Save all ORA results from a nested run
#   save_all_gsea_results()    - Save all GSEA results from a run
#   build_pathway_driver_table() - Flat per-pathway table w/ driver genes (Panther-style)
#   write_pathway_driver_workbook() - Excel workbook of driver tables per contrast
#
# Dependencies (install in your conda env):
#   clusterProfiler, enrichplot, msigdbr, ggplot2, dplyr, openxlsx,
#   tidyr, ggrepel, RColorBrewer
#
# Optional dependencies for visualization:
#   ggupset (for upsetplot), patchwork (for combining plots)
# ============================================================================


# --------------------------------------------------------------------------
# load_geneset_collection()
# --------------------------------------------------------------------------
#' Load one MSigDB gene set collection in TERM2GENE format
#'
#' @param collection MSigDB collection ID. Common ones:
#'   "H"  - Hallmark gene sets (50 broad themes)
#'   "C2" - Curated (incl. KEGG, Reactome subcollections)
#'   "C5" - Ontology (GO subcollections)
#'   "C7" - Immunologic signatures
#' @param subcollection Optional subcollection. Examples:
#'   For C2: "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:REACTOME", "CP:WIKIPATHWAYS"
#'   For C5: "GO:BP", "GO:MF", "GO:CC"
#' @param species Species name as msigdbr expects ("Mus musculus", "Homo sapiens", etc.)
#' @return data.frame with columns gs_name, gene_symbol (TERM2GENE format)
load_geneset_collection <- function(collection, subcollection = NULL,
                                    species = "Mus musculus") {
  args <- list(species = species, collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection

  df <- do.call(msigdbr::msigdbr, args)
  if (nrow(df) == 0) {
    stop("No gene sets returned for collection=", collection,
         ", subcollection=", subcollection, ", species=", species)
  }

  out <- df %>% dplyr::select(gs_name, gene_symbol)
  attr(out, "collection")    <- collection
  attr(out, "subcollection") <- subcollection
  attr(out, "label")         <- if (is.null(subcollection)) collection
                                else paste0(collection, "_", subcollection)
  out
}


# --------------------------------------------------------------------------
# load_all_collections()
# --------------------------------------------------------------------------
#' Load multiple gene set collections at once
#'
#' @param specs A named list. Each element is a list with `collection` and
#'   optional `subcollection`. The names of the list become the labels.
#'   Default loads our "full stack": GO:BP, Hallmark, KEGG, Reactome.
#' @param species Species name. Default "Mus musculus".
#' @return Named list of TERM2GENE data.frames
load_all_collections <- function(specs = NULL, species = "Mus musculus") {
  if (is.null(specs)) {
    specs <- list(
      Hallmark = list(collection = "H"),
      GO_BP    = list(collection = "C5", subcollection = "GO:BP"),
      KEGG     = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
      Reactome = list(collection = "C2", subcollection = "CP:REACTOME")
    )
  }

  cat("Loading", length(specs), "gene set collections for", species, "...\n")
  out <- lapply(names(specs), function(label) {
    spec <- specs[[label]]
    cat("  ", label, "...")
    coll <- load_geneset_collection(
      collection = spec$collection,
      subcollection = spec$subcollection,
      species = species
    )
    cat(" loaded", length(unique(coll$gs_name)), "gene sets\n")
    coll
  })
  names(out) <- names(specs)
  out
}


# --------------------------------------------------------------------------
# prepare_ora_input()
# --------------------------------------------------------------------------
#' Build ORA-ready gene lists from a DESeq2 results table
#'
#' Returns three foreground gene lists (UP, DOWN, ALL) plus the background.
#'
#' @param results_df A data.frame with at minimum: gene_name, log2FoldChange, padj.
#'   Typically loaded from DESeq2 *_full_results.csv outputs.
#' @param padj_threshold padj cutoff. Default 0.05.
#' @param lfc_threshold |log2FC| cutoff. Default 1 (2-fold). Use 0.585 for 1.5-fold.
#' @param symbol_col Column name for gene symbols. Default "gene_name".
#' @return List with: up, down, all (foreground gene vectors), background, and counts
prepare_ora_input <- function(results_df,
                              padj_threshold = 0.05,
                              lfc_threshold = 1,
                              symbol_col = "gene_name") {

  if (!symbol_col %in% colnames(results_df)) {
    stop("Symbol column '", symbol_col, "' not found in results_df.")
  }

  # Background = all genes with non-NA padj (i.e., genes that were tested)
  bg <- results_df %>%
    dplyr::filter(!is.na(padj),
                  !is.na(.data[[symbol_col]]),
                  .data[[symbol_col]] != "") %>%
    dplyr::pull(!!rlang::sym(symbol_col)) %>%
    unique()

  # Significant set
  sig <- results_df %>%
    dplyr::filter(!is.na(padj),
                  padj < padj_threshold,
                  abs(log2FoldChange) >= lfc_threshold,
                  !is.na(.data[[symbol_col]]),
                  .data[[symbol_col]] != "")

  up   <- sig %>% dplyr::filter(log2FoldChange > 0) %>%
                  dplyr::pull(!!rlang::sym(symbol_col)) %>% unique()
  down <- sig %>% dplyr::filter(log2FoldChange < 0) %>%
                  dplyr::pull(!!rlang::sym(symbol_col)) %>% unique()
  all  <- unique(c(up, down))

  cat("ORA input prepared:\n")
  cat("  Background genes:    ", length(bg), "\n")
  cat("  Significant (any):   ", length(all), "\n")
  cat("  Significant (UP):    ", length(up), "\n")
  cat("  Significant (DOWN):  ", length(down), "\n")

  list(up = up, down = down, all = all, background = bg,
       n_bg = length(bg), n_up = length(up),
       n_down = length(down), n_all = length(all))
}


# --------------------------------------------------------------------------
# prepare_gsea_input()
# --------------------------------------------------------------------------
#' Build a ranked named vector for GSEA from a DESeq2 results table
#'
#' Default ranking: log2FoldChange (sorted high to low). Alternative metrics
#' available via `metric` argument.
#'
#' @param results_df Same as prepare_ora_input().
#' @param metric "lfc" (default), "stat" (-log10(p) * sign(LFC)), or "shrunk_lfc"
#'   (same as "lfc" but emphasizes named choice).
#' @param symbol_col Column name for gene symbols. Default "gene_name".
#' @return Named numeric vector, sorted descending. Names = gene symbols.
prepare_gsea_input <- function(results_df,
                               metric = "lfc",
                               symbol_col = "gene_name") {

  results_df <- results_df %>%
    dplyr::filter(!is.na(log2FoldChange),
                  !is.na(.data[[symbol_col]]),
                  .data[[symbol_col]] != "")

  if (metric == "lfc" || metric == "shrunk_lfc") {
    results_df$rank_metric <- results_df$log2FoldChange
  } else if (metric == "stat") {
    if (!"pvalue" %in% colnames(results_df)) {
      stop("'stat' metric requires a 'pvalue' column.")
    }
    results_df <- results_df %>% dplyr::filter(!is.na(pvalue))
    results_df$rank_metric <- -log10(results_df$pvalue) *
                              sign(results_df$log2FoldChange)
  } else {
    stop("Unknown metric: ", metric, ". Use 'lfc' or 'stat'.")
  }

  # For duplicate gene symbols, keep the one with the largest absolute rank
  out <- results_df %>%
    dplyr::group_by(.data[[symbol_col]]) %>%
    dplyr::slice_max(abs(rank_metric), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rank_metric)) %>%
    dplyr::select(!!rlang::sym(symbol_col), rank_metric) %>%
    tibble::deframe()

  cat("GSEA input prepared:\n")
  cat("  Total ranked genes:", length(out), "\n")
  cat("  Range:            ", round(min(out), 2), "to", round(max(out), 2), "\n")

  out
}


# --------------------------------------------------------------------------
# run_ora()
# --------------------------------------------------------------------------
#' Run a single ORA analysis with our standard settings
#'
#' @param genes Foreground gene vector (gene symbols).
#' @param TERM2GENE TERM2GENE data.frame from load_geneset_collection().
#' @param universe Background gene vector. Critical for unbiased results.
#' @param padj_cutoff BH-adjusted p-value cutoff. Default 0.05.
#' @param min_gs_size Minimum overlap between gene set and background. Default 10.
#' @param max_gs_size Maximum gene set size. Default 500 (filters out massive sets).
#' @return enrichResult object (clusterProfiler S4 class). Returns NULL if no genes.
run_ora <- function(genes, TERM2GENE, universe,
                    padj_cutoff = 0.05,
                    min_gs_size = 10,
                    max_gs_size = 500) {

  if (length(genes) == 0) {
    cat("  (skipped - no genes in foreground)\n")
    return(NULL)
  }

  res <- tryCatch(
    clusterProfiler::enricher(
      gene          = genes,
      universe      = universe,
      TERM2GENE     = TERM2GENE,
      pvalueCutoff  = 1,    # we'll filter on padj after
      qvalueCutoff  = 1,
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      cat("  ORA failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat("  (no significant pathways)\n")
    return(res)
  }

  # Filter the @result slot to padj_cutoff while preserving the object structure
  res@result <- res@result[res@result$p.adjust < padj_cutoff, , drop = FALSE]
  cat("  ", nrow(as.data.frame(res)), "significant pathways at padj <", padj_cutoff, "\n")
  res
}


# --------------------------------------------------------------------------
# run_ora_all_collections()
# --------------------------------------------------------------------------
#' Run ORA across multiple gene set collections AND directions in one call
#'
#' Produces a nested list: results[[collection]][[direction]] = enrichResult
#' Directions = "all", "up", "down".
#'
#' @param ora_input Output of prepare_ora_input().
#' @param collections Named list of TERM2GENE data.frames from load_all_collections().
#' @param directions Which directions to test. Default c("all", "up", "down").
#' @param ... Passed to run_ora().
#' @return Nested list: results[[collection_name]][[direction]] = enrichResult
run_ora_all_collections <- function(ora_input, collections,
                                    directions = c("all", "up", "down"),
                                    ...) {

  results <- list()
  for (coll_name in names(collections)) {
    cat("\n=== Collection:", coll_name, "===\n")
    coll_results <- list()
    for (dir in directions) {
      genes <- ora_input[[dir]]
      cat("  Direction:", dir, "(", length(genes), "genes) ... ")
      coll_results[[dir]] <- run_ora(
        genes = genes,
        TERM2GENE = collections[[coll_name]],
        universe = ora_input$background,
        ...
      )
    }
    results[[coll_name]] <- coll_results
  }
  results
}


# --------------------------------------------------------------------------
# run_gsea()
# --------------------------------------------------------------------------
#' Run a single GSEA analysis with our standard settings
#'
#' @param gene_list Named numeric vector from prepare_gsea_input().
#' @param TERM2GENE TERM2GENE data.frame.
#' @param padj_cutoff Cutoff (default 0.05); only filters reported results, not the test.
#' @param min_gs_size Minimum gene set size. Default 10.
#' @param max_gs_size Maximum gene set size. Default 500.
#' @param seed Random seed for permutations. Default 42.
#' @return gseaResult object, or NULL if no significant results.
run_gsea <- function(gene_list, TERM2GENE,
                     padj_cutoff = 0.05,
                     min_gs_size = 10,
                     max_gs_size = 500,
                     seed = 42) {

  res <- tryCatch(
    clusterProfiler::GSEA(
      geneList      = gene_list,
      TERM2GENE     = TERM2GENE,
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      pvalueCutoff  = 1,    # we'll filter after
      pAdjustMethod = "BH",
      seed          = seed,
      nPermSimple   = 10000,
      verbose       = FALSE
    ),
    error = function(e) {
      cat("  GSEA failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat("  (no significant gene sets)\n")
    return(res)
  }

  res@result <- res@result[res@result$p.adjust < padj_cutoff, , drop = FALSE]
  cat("  ", nrow(as.data.frame(res)), "significant gene sets at padj <", padj_cutoff, "\n")
  res
}


# --------------------------------------------------------------------------
# run_gsea_all_collections()
# --------------------------------------------------------------------------
#' Run GSEA across multiple gene set collections in one call
#'
#' @param gene_list Named numeric vector from prepare_gsea_input().
#' @param collections Named list from load_all_collections().
#' @param ... Passed to run_gsea().
#' @return Named list: results[[collection_name]] = gseaResult
run_gsea_all_collections <- function(gene_list, collections, ...) {
  results <- list()
  for (coll_name in names(collections)) {
    cat("\n=== Collection:", coll_name, "===\n")
    cat("  Running GSEA ... ")
    results[[coll_name]] <- run_gsea(
      gene_list = gene_list,
      TERM2GENE = collections[[coll_name]],
      ...
    )
  }
  results
}


# --------------------------------------------------------------------------
# plot_ora_dotplot()
# --------------------------------------------------------------------------
#' Standard dotplot of top ORA results
#'
#' @param ora_result enrichResult object from run_ora().
#' @param show_n Number of top pathways to show. Default 20.
#' @param title Plot title.
#' @return ggplot object, or NULL if no results.
plot_ora_dotplot <- function(ora_result, show_n = 20, title = NULL) {
  if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
    return(NULL)
  }
  p <- enrichplot::dotplot(ora_result,
                            showCategory = show_n,
                            x = "GeneRatio",
                            color = "p.adjust",
                            size = "Count")
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}


# --------------------------------------------------------------------------
# plot_ora_emap()
# --------------------------------------------------------------------------
#' Enrichment map: pathways clustered by gene overlap
#'
#' @param ora_result enrichResult.
#' @param show_n Number of pathways. Default 30.
#' @param title Plot title.
#' @return ggplot object, or NULL if too few results.
plot_ora_emap <- function(ora_result, show_n = 30, title = NULL) {
  if (is.null(ora_result) || nrow(as.data.frame(ora_result)) < 2) {
    return(NULL)
  }
  sim <- enrichplot::pairwise_termsim(ora_result)
  p <- enrichplot::emapplot(sim, showCategory = show_n)
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}


# --------------------------------------------------------------------------
# plot_gsea_dotplot()
# --------------------------------------------------------------------------
#' Dotplot of top GSEA results, separately for activated and suppressed pathways
#'
#' @param gsea_result gseaResult object.
#' @param show_n Number of pathways per direction.
#' @param title Plot title.
#' @return ggplot object.
plot_gsea_dotplot <- function(gsea_result, show_n = 15, title = NULL) {
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    return(NULL)
  }
  p <- tryCatch(
    enrichplot::dotplot(gsea_result, showCategory = show_n,
                         split = ".sign", x = "NES",
                         color = "p.adjust") +
      ggplot2::facet_grid(. ~ .sign),
    error = function(e) {
      enrichplot::dotplot(gsea_result, showCategory = show_n, x = "NES",
                           color = "p.adjust")
    }
  )
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}


# --------------------------------------------------------------------------
# plot_gsea_top_n()
# --------------------------------------------------------------------------
#' Print the top N GSEA enrichment plots one after another
#'
#' @param gsea_result gseaResult object.
#' @param n Number of top gene sets to plot.
#' @param save_dir Optional directory to save individual PDFs.
#' @param prefix Filename prefix if saving.
plot_gsea_top_n <- function(gsea_result, n = 5, save_dir = NULL, prefix = "gsea") {
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    cat("(no results to plot)\n")
    return(invisible(NULL))
  }
  top_n <- as.data.frame(gsea_result) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n)

  for (i in seq_len(nrow(top_n))) {
    p <- enrichplot::gseaplot2(
      gsea_result,
      geneSetID = top_n$ID[i],
      title = paste0("#", i, " | ", top_n$Description[i],
                     "\nNES = ", round(top_n$NES[i], 3),
                     " | FDR = ", signif(top_n$p.adjust[i], 3)),
      pvalue_table = TRUE
    )
    print(p)

    if (!is.null(save_dir)) {
      ggplot2::ggsave(
        file.path(save_dir, paste0(prefix, "_top", i, "_",
                                    gsub("[^A-Za-z0-9]+", "_",
                                         substr(top_n$Description[i], 1, 50)),
                                    ".pdf")),
        plot = p, width = 10, height = 8
      )
    }
  }
}


# --------------------------------------------------------------------------
# compare_ora_dotplot()
# --------------------------------------------------------------------------
#' Compare ORA results across multiple contrasts side-by-side as a dotplot
#'
#' Takes a named list of enrichResult objects, picks the union of top N pathways
#' from each, and shows them all as a dotplot with contrasts on the X axis.
#'
#' @param ora_list Named list of enrichResult objects (names become column labels).
#' @param top_n_per_contrast How many top pathways from each contrast to include.
#' @param title Plot title.
#' @return ggplot object.
compare_ora_dotplot <- function(ora_list, top_n_per_contrast = 10, title = NULL) {

  # Collect top N pathways from each contrast
  top_pathways <- character(0)
  combined <- list()
  for (nm in names(ora_list)) {
    res <- ora_list[[nm]]
    if (is.null(res) || nrow(as.data.frame(res)) == 0) next
    df <- as.data.frame(res) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n_per_contrast)
    top_pathways <- union(top_pathways, df$Description)
    combined[[nm]] <- as.data.frame(res) %>%
      dplyr::mutate(contrast = nm)
  }

  if (length(combined) == 0 || length(top_pathways) == 0) {
    cat("No pathways to plot.\n")
    return(NULL)
  }

  full_df <- dplyr::bind_rows(combined) %>%
    dplyr::filter(Description %in% top_pathways) %>%
    dplyr::mutate(
      neglog10padj = -log10(p.adjust),
      gene_ratio_num = sapply(GeneRatio, function(s) {
        parts <- as.numeric(strsplit(s, "/")[[1]])
        if (length(parts) == 2) parts[1] / parts[2] else NA
      }),
      contrast = factor(contrast, levels = names(ora_list))
    )

  # Truncate long pathway names
  full_df$Description_short <- ifelse(nchar(full_df$Description) > 60,
                                       paste0(substr(full_df$Description, 1, 57), "..."),
                                       full_df$Description)

  # Order pathways by max significance across contrasts
  pathway_order <- full_df %>%
    dplyr::group_by(Description_short) %>%
    dplyr::summarise(max_sig = max(neglog10padj, na.rm = TRUE)) %>%
    dplyr::arrange(max_sig) %>%
    dplyr::pull(Description_short)
  full_df$Description_short <- factor(full_df$Description_short, levels = pathway_order)

  p <- ggplot2::ggplot(full_df,
                        ggplot2::aes(x = contrast, y = Description_short,
                                      color = neglog10padj, size = gene_ratio_num)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "blue", high = "red",
                                   name = "-log10(padj)") +
    ggplot2::scale_size_continuous(name = "GeneRatio") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = title)
  p
}


# --------------------------------------------------------------------------
# save_enrichment_results()
# --------------------------------------------------------------------------
#' Save ORA or GSEA results to CSV + Excel
#'
#' @param result_obj enrichResult or gseaResult, or NULL.
#' @param output_dir Output directory.
#' @param filename_prefix Prefix for filenames (no extension).
#' @return Invisibly the data.frame written, or NULL if nothing.
save_enrichment_results <- function(result_obj, output_dir, filename_prefix) {
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {
    cat("  (no results to save for ", filename_prefix, ")\n", sep = "")
    return(invisible(NULL))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  df <- as.data.frame(result_obj)
  df$minuslog10padj <- -log10(df$p.adjust)

  csv_path <- file.path(output_dir, paste0(filename_prefix, ".csv"))
  xlsx_path <- file.path(output_dir, paste0(filename_prefix, ".xlsx"))

  write.csv(df, csv_path, row.names = FALSE)
  openxlsx::write.xlsx(df, xlsx_path, asTable = TRUE,
                        tableStyle = "TableStyleMedium2", firstRow = TRUE)

  cat("  Saved", filename_prefix, "(", nrow(df), "rows)\n")
  invisible(df)
}


# --------------------------------------------------------------------------
# save_all_ora_results()
# --------------------------------------------------------------------------
#' Save all ORA results from a nested run_ora_all_collections() output
#'
#' @param ora_results Output of run_ora_all_collections().
#' @param output_dir Where to save files.
#' @param contrast_name Used as filename prefix.
save_all_ora_results <- function(ora_results, output_dir, contrast_name) {
  for (coll_name in names(ora_results)) {
    for (dir in names(ora_results[[coll_name]])) {
      save_enrichment_results(
        ora_results[[coll_name]][[dir]],
        output_dir = output_dir,
        filename_prefix = paste0(contrast_name, "_ORA_", coll_name, "_", dir)
      )
    }
  }
}

#' Save all GSEA results from a run_gsea_all_collections() output
save_all_gsea_results <- function(gsea_results, output_dir, contrast_name) {
  for (coll_name in names(gsea_results)) {
    save_enrichment_results(
      gsea_results[[coll_name]],
      output_dir = output_dir,
      filename_prefix = paste0(contrast_name, "_GSEA_", coll_name)
    )
  }
}


# --------------------------------------------------------------------------
# build_pathway_driver_table()
# --------------------------------------------------------------------------
#' Build a flat per-pathway table with driver genes for one enrichment result
#'
#' Mirrors the supplemental-table format used in many RNA-seq papers
#' (e.g. the Panther GO Biological Processes tables): one row per pathway,
#' with pathway metadata and the actual gene symbols that drove the
#' enrichment listed as a comma-separated string.
#'
#' Handles both GSEA and ORA results:
#'   - GSEA: driver genes are taken from `core_enrichment` (leading-edge)
#'           direction inferred from sign of NES
#'   - ORA:  driver genes are taken from `geneID` (overlap genes)
#'           direction taken from the `direction` argument
#'
#' @param result_obj enrichResult or gseaResult object from clusterProfiler.
#'   Pass NULL or an empty result and you'll get back an empty data frame.
#' @param method "GSEA" or "ORA". Determines which fields to read.
#' @param direction For ORA only: "up", "down", or "all". Adds a Direction column.
#'   Ignored for GSEA (derived from NES sign automatically).
#' @param collection_name Optional collection name to include as a column.
#' @param max_pathways Cap the number of pathways included. Default 50.
#'   Sorted by p.adjust ascending before capping.
#' @return A data.frame with columns:
#'   Direction, Collection, ID, Description, NES (GSEA) or FoldEnrichment (ORA),
#'   p.adjust, minuslog10padj, GeneCount, Genes
#'   For empty inputs, returns a data frame with zero rows.
build_pathway_driver_table <- function(result_obj,
                                       method = c("GSEA", "ORA"),
                                       direction = NULL,
                                       collection_name = NA_character_,
                                       max_pathways = 50) {
  method <- match.arg(method)

  # Handle empty/null inputs gracefully
  if (is.null(result_obj) || nrow(as.data.frame(result_obj)) == 0) {
    return(data.frame(
      Direction = character(0),
      Collection = character(0),
      ID = character(0),
      Description = character(0),
      NES_or_FoldEnrichment = numeric(0),
      p.adjust = numeric(0),
      minuslog10padj = numeric(0),
      GeneCount = integer(0),
      Genes = character(0),
      stringsAsFactors = FALSE
    ))
  }

  df <- as.data.frame(result_obj)
  df <- df[order(df$p.adjust), , drop = FALSE]
  if (nrow(df) > max_pathways) {
    df <- df[seq_len(max_pathways), , drop = FALSE]
  }

  # Method-specific columns
  if (method == "GSEA") {
    # Direction from NES sign
    dir_vec <- ifelse(df$NES > 0, "positive", "negative")
    effect_size <- df$NES
    effect_name <- "NES"
    # Leading-edge genes
    gene_field <- df$core_enrichment
  } else {  # ORA
    if (is.null(direction)) {
      stop("`direction` argument is required for ORA results (\"up\", \"down\", or \"all\").")
    }
    dir_vec <- direction
    # Compute fold enrichment from GeneRatio / BgRatio (Panther-style metric)
    parse_ratio <- function(x) {
      parts <- strsplit(x, "/", fixed = TRUE)
      sapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]))
    }
    fold_enr <- parse_ratio(df$GeneRatio) / parse_ratio(df$BgRatio)
    effect_size <- fold_enr
    effect_name <- "FoldEnrichment"
    gene_field <- df$geneID
  }

  # Convert slash-separated gene strings to comma-separated (Panther style)
  genes_pretty <- vapply(gene_field, function(g) {
    if (is.na(g) || g == "") return("")
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    paste(parts, collapse = ", ")
  }, character(1))
  gene_counts <- vapply(gene_field, function(g) {
    if (is.na(g) || g == "") return(0L)
    length(strsplit(g, "/", fixed = TRUE)[[1]])
  }, integer(1))

  out <- data.frame(
    Direction      = dir_vec,
    Collection     = rep(collection_name, nrow(df)),
    ID             = df$ID,
    Description    = df$Description,
    NES_or_FoldEnrichment = effect_size,
    p.adjust       = df$p.adjust,
    minuslog10padj = -log10(df$p.adjust),
    GeneCount      = gene_counts,
    Genes          = genes_pretty,
    stringsAsFactors = FALSE
  )
  # Rename the effect-size column to its actual meaning for clarity
  colnames(out)[colnames(out) == "NES_or_FoldEnrichment"] <- effect_name
  rownames(out) <- NULL

  out
}


# --------------------------------------------------------------------------
# write_pathway_driver_workbook()
# --------------------------------------------------------------------------
#' Write a Panther-style Excel workbook of pathway-driver gene tables
#'
#' Produces ONE Excel workbook per contrast, with separate sheets for each
#' collection × method combination. Each sheet is a flat table mirroring
#' the supplemental-table format from many RNA-seq papers.
#'
#' Sheets are named like "GSEA_Hallmark", "ORA_GOBP_up", "ORA_GOBP_down".
#'
#' @param gsea_results Optional. Output of run_gsea_all_collections().
#'   Named list: results[[collection]] = gseaResult.
#' @param ora_results Optional. Output of run_ora_all_collections().
#'   Nested: results[[collection]][[direction]] = enrichResult.
#' @param output_path Path to output .xlsx file.
#' @param max_pathways_per_sheet Cap on rows per sheet. Default 50.
#' @return Invisibly, the output path.
write_pathway_driver_workbook <- function(gsea_results = NULL,
                                          ora_results  = NULL,
                                          output_path,
                                          max_pathways_per_sheet = 50) {

  if (is.null(gsea_results) && is.null(ora_results)) {
    stop("Provide at least one of gsea_results or ora_results.")
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  wb <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(textDecoration = "bold",
                                         fgFill = "#D9E1F2",
                                         border = "bottom")
  any_sheets <- FALSE

  # ---- GSEA sheets ----
  if (!is.null(gsea_results)) {
    for (coll_name in names(gsea_results)) {
      tbl <- build_pathway_driver_table(
        gsea_results[[coll_name]],
        method = "GSEA",
        collection_name = coll_name,
        max_pathways = max_pathways_per_sheet
      )
      if (nrow(tbl) == 0) {
        cat("  GSEA / ", coll_name, ": no pathways to write\n", sep = "")
        next
      }
      sheet_name <- paste0("GSEA_", coll_name)
      # Excel sheet names cap at 31 chars
      sheet_name <- substr(gsub("[/\\\\?*\\[\\]:]", "_", sheet_name), 1, 31)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, tbl, headerStyle = header_style)
      openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)),
                              widths = "auto")
      # Wrap the Genes column (last column)
      openxlsx::addStyle(wb, sheet_name,
                         style = openxlsx::createStyle(wrapText = TRUE,
                                                       valign = "top"),
                         rows = 2:(nrow(tbl) + 1),
                         cols = ncol(tbl), gridExpand = TRUE)
      openxlsx::setColWidths(wb, sheet_name, cols = ncol(tbl), widths = 80)
      openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
      cat("  Wrote sheet:", sheet_name, "(", nrow(tbl), "pathways)\n")
      any_sheets <- TRUE
    }
  }

  # ---- ORA sheets ----
  if (!is.null(ora_results)) {
    for (coll_name in names(ora_results)) {
      for (dir_name in names(ora_results[[coll_name]])) {
        tbl <- build_pathway_driver_table(
          ora_results[[coll_name]][[dir_name]],
          method = "ORA",
          direction = dir_name,
          collection_name = coll_name,
          max_pathways = max_pathways_per_sheet
        )
        if (nrow(tbl) == 0) {
          cat("  ORA / ", coll_name, " / ", dir_name,
              ": no pathways to write\n", sep = "")
          next
        }
        sheet_name <- paste0("ORA_", coll_name, "_", dir_name)
        sheet_name <- substr(gsub("[/\\\\?*\\[\\]:]", "_", sheet_name), 1, 31)
        openxlsx::addWorksheet(wb, sheet_name)
        openxlsx::writeData(wb, sheet_name, tbl, headerStyle = header_style)
        openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(tbl)),
                                widths = "auto")
        openxlsx::addStyle(wb, sheet_name,
                           style = openxlsx::createStyle(wrapText = TRUE,
                                                         valign = "top"),
                           rows = 2:(nrow(tbl) + 1),
                           cols = ncol(tbl), gridExpand = TRUE)
        openxlsx::setColWidths(wb, sheet_name, cols = ncol(tbl), widths = 80)
        openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
        cat("  Wrote sheet:", sheet_name, "(", nrow(tbl), "pathways)\n")
        any_sheets <- TRUE
      }
    }
  }

  if (!any_sheets) {
    cat("  WARNING: no pathways found in any input - skipping workbook write\n")
    return(invisible(NULL))
  }

  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  cat("\nSaved pathway-driver workbook:", output_path, "\n")
  invisible(output_path)
}


# --------------------------------------------------------------------------
# plot_gsea_horizontal_bars()
# --------------------------------------------------------------------------
#' Publication-style horizontal bar plot of GSEA pathways (de Man Aging Cell style)
#'
#' Bar length = signed -log10(padj), where the sign follows the NES sign:
#' positive NES = bar extending right (up-regulated pathway), negative NES =
#' bar extending left (down-regulated pathway). NES is shown as a text
#' annotation at the end of each bar, matching the "Fold Enrichment" numbers
#' in the de Man figure. Color is direction (red = up, blue = down).
#'
#' @param gsea_result A gseaResult object from clusterProfiler::GSEA(), or NULL.
#' @param n_up Number of top up-regulated (positive NES) pathways to show. Default 5.
#' @param n_down Number of top down-regulated (negative NES) pathways to show. Default 5.
#' @param padj_cutoff Only include pathways at or below this padj. Default 0.05.
#' @param title Plot title.
#' @param wrap_width Wrap pathway names at this character width for readability. Default 40.
#' @param show_nes Logical; whether to print NES numbers at the end of each bar. Default TRUE.
#' @param up_color Hex color for positive-NES bars. Default brick red.
#' @param down_color Hex color for negative-NES bars. Default steel blue.
#' @return A ggplot object, or NULL if no significant pathways.
plot_gsea_horizontal_bars <- function(gsea_result,
                                      n_up = 5, n_down = 5,
                                      padj_cutoff = 0.05,
                                      title = NULL,
                                      wrap_width = 40,
                                      show_nes = TRUE,
                                      up_color = "#C0392B",
                                      down_color = "#2874A6") {

  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    message("No GSEA results to plot.")
    return(NULL)
  }

  df <- as.data.frame(gsea_result)
  df <- df[!is.na(df$p.adjust) & df$p.adjust <= padj_cutoff, , drop = FALSE]
  if (nrow(df) == 0) {
    message("No pathways pass padj_cutoff = ", padj_cutoff, ".")
    return(NULL)
  }

  # Pick top n_up positive-NES and top n_down negative-NES, sorted by padj
  up   <- df[df$NES > 0, , drop = FALSE]
  down <- df[df$NES < 0, , drop = FALSE]
  up   <- up[order(up$p.adjust), , drop = FALSE]
  down <- down[order(down$p.adjust), , drop = FALSE]
  if (nrow(up)   > n_up)   up   <- up[seq_len(n_up), , drop = FALSE]
  if (nrow(down) > n_down) down <- down[seq_len(n_down), , drop = FALSE]
  sel <- rbind(up, down)

  if (nrow(sel) == 0) {
    message("No pathways to plot after filtering by direction.")
    return(NULL)
  }

  # Signed -log10(padj): sign comes from NES direction
  sel$signed_logp <- sign(sel$NES) * -log10(sel$p.adjust)
  sel$direction <- ifelse(sel$NES > 0, "Up", "Down")

  # Clean up MSigDB pathway names: drop "HALLMARK_" / "GOBP_" / "REACTOME_"
  # prefixes and replace underscores with spaces for readability
  sel$pretty_name <- gsub("^[A-Z]+_", "", sel$Description)
  sel$pretty_name <- gsub("_", " ", sel$pretty_name)
  sel$pretty_name <- tools::toTitleCase(tolower(sel$pretty_name))
  # Wrap long names
  sel$pretty_name <- vapply(sel$pretty_name, function(x) {
    paste(strwrap(x, width = wrap_width), collapse = "\n")
  }, character(1))

  # Order y-axis: most positive at top, most negative at bottom
  sel <- sel[order(sel$signed_logp), , drop = FALSE]
  sel$pretty_name <- factor(sel$pretty_name, levels = sel$pretty_name)

  # NES annotation, with alignment depending on bar direction:
  # positive bars: annotation placed just past the right end (hjust = -0.15)
  # negative bars: annotation placed just past the left end  (hjust = 1.15)
  sel$nes_label <- sprintf("%.2f", sel$NES)
  sel$nes_hjust <- ifelse(sel$signed_logp > 0, -0.15, 1.15)

  # Build plot
  x_range <- range(sel$signed_logp, na.rm = TRUE)
  # Pad so the NES labels don't get clipped
  x_pad <- diff(x_range) * 0.18
  x_lim <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  p <- ggplot(sel, aes(x = signed_logp, y = pretty_name, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c("Up" = up_color, "Down" = down_color),
                      breaks = c("Up", "Down")) +
    scale_x_continuous(limits = x_lim,
                       labels = function(x) abs(x)) +
    labs(
      x = expression(-log[10] * "(p.adjust)  ("*phantom() %<-% "down  /  up" %->% phantom()*")"),
      y = NULL,
      fill = "Direction",
      title = title
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title = element_text(face = "bold")
    )

  if (show_nes) {
    p <- p + geom_text(aes(label = nes_label, hjust = nes_hjust),
                       size = 3.4, color = "black")
  }

  p
}
