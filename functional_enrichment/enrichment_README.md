# Functional Enrichment Helpers

Reusable helper functions for ORA + GSEA analysis of any DESeq2 results,
using clusterProfiler and MSigDB gene sets.

## Files

```
functional_enrichment/
├── README.md            <- this file
└── enrichment_utils.R   <- helper functions
```

## How to use

In any analysis notebook, source the file with the absolute path:

```r
source("~/project_pi_epm38/jfd42/pipelines/functional_enrichment/enrichment_utils.R")
```

Then use the helpers as documented in each function's header comment.

## Function inventory

| Category | Function | Purpose |
|----------|----------|---------|
| Loading | `load_geneset_collection()` | Load one MSigDB collection |
| Loading | `load_all_collections()` | Load multiple at once (default = full stack) |
| Input prep | `prepare_ora_input()` | Build UP/DOWN/ALL gene lists + background from DESeq2 results |
| Input prep | `prepare_gsea_input()` | Build ranked named vector for GSEA |
| Run | `run_ora()` | Single ORA test |
| Run | `run_ora_all_collections()` | All collections × directions |
| Run | `run_gsea()` | Single GSEA test |
| Run | `run_gsea_all_collections()` | All collections |
| Plot | `plot_ora_dotplot()` | Standard dotplot |
| Plot | `plot_ora_emap()` | Enrichment map |
| Plot | `plot_gsea_dotplot()` | GSEA dotplot |
| Plot | `plot_gsea_top_n()` | Grid of top N enrichment plots |
| Plot | `compare_ora_dotplot()` | Side-by-side multi-contrast comparison |
| Save | `save_enrichment_results()` | Single CSV + Excel writeout |
| Save | `save_all_ora_results()` | Bulk save from nested list |
| Save | `save_all_gsea_results()` | Bulk save GSEA results |

## Default gene set collections

`load_all_collections()` with no arguments loads:

- **Hallmark** (50 broad themes, fast big-picture overview)
- **GO:BP** (Gene Ontology Biological Process — fine-grained, the most informative for biology)
- **KEGG** (curated pathway database — KEGG_LEGACY subcollection)
- **Reactome** (curated pathways, complementary to KEGG)

Override by passing your own `specs` list. See the function header for format.

## Species

Default is `Mus musculus`. Override with the `species` argument. To see what
species msigdbr supports:

```r
msigdbr::msigdbr_species()
```

For organisms not in msigdbr, you can convert from human/mouse via orthologs
using biomaRt or gprofiler2, then pass a custom TERM2GENE data.frame to
`run_ora()` directly.

## Dependencies

Should already be in your bulkrnaseq conda environment:
- clusterProfiler
- enrichplot
- msigdbr
- ggplot2 / dplyr / tidyr (via tidyverse)
- openxlsx
- ggrepel

Optional:
- ggupset (for upsetplot)
- patchwork (for combining plots)
