# ============================================================
# Project 2: Alzheimer's Neurobioinformatics
# Stage 5: Signed Wald rankings + GO BP and KEGG GSEA
# ============================================================

# Install required Bioconductor packages only if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  BiocManager::install(
    missing_packages,
    ask = FALSE,
    update = FALSE
  )
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(20261021)

# Define directories cleanly as strings
output_dir <- "results/gsea"
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Convert DESeq2 results into a ranked Entrez-ID list
make_ranked_list <- function(res_obj) {

  rank_df <- data.frame(
    SYMBOL = rownames(res_obj),
    stat = res_obj$stat,
    stringsAsFactors = FALSE
  )

  rank_df <- rank_df[
    is.finite(rank_df$stat) &
      !is.na(rank_df$SYMBOL) &
      nzchar(rank_df$SYMBOL),
    ,
    drop = FALSE
  ]

  # Map gene symbols to Entrez identifiers
  id_map <- clusterProfiler::bitr(
    unique(rank_df$SYMBOL),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  mapped <- merge(
    rank_df,
    id_map,
    by = "SYMBOL"
  )

  # Keep one entry per Entrez ID.
  mapped <- mapped[
    order(abs(mapped$stat), decreasing = TRUE),
    ,
    drop = FALSE
  ]

  mapped <- mapped[
    !duplicated(mapped$ENTREZID),
    ,
    drop = FALSE
  ]

  gene_list <- mapped$stat
  names(gene_list) <- mapped$ENTREZID

  sort(gene_list, decreasing = TRUE)
}

rank_AD_Old <- make_ranked_list(res_AD_Old)
rank_Old_Young <- make_ranked_list(res_Old_Young)

cat("Mapped genes, AD vs Old:", length(rank_AD_Old), "\n")
cat("Mapped genes, Old vs Young:", length(rank_Old_Young), "\n")

# FIXED: Wrapped directory targets in quotes/pre-defined variables
saveRDS(
  rank_AD_Old,
  file.path(output_dir, "ranked_genes_AD_vs_Old.rds")
)

saveRDS(
  rank_Old_Young,
  file.path(output_dir, "ranked_genes_Old_vs_Young.rds")
)

# 3. Run GO Biological Process and KEGG GSEA
run_gsea <- function(gene_list, label) {

  cat("\nRunning GSEA for:", label, "\n")

  go_bp <- clusterProfiler::gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    verbose = FALSE
  )

  kegg <- tryCatch(
    clusterProfiler::gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      keyType = "ncbi-geneid",
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      verbose = FALSE
    ),
    error = function(e) {
      message("KEGG GSEA could not complete: ", conditionMessage(e))
      NULL
    }
  )

  # Save and display pathway results
  save_and_show <- function(result, database) {

    if (is.null(result)) {
      return(invisible(NULL))
    }

    result_df <- as.data.frame(result)

    out_file <- file.path(
      output_dir,
      paste0(database, "_", label, ".csv")
    )

    write.csv(
      result_df,
      out_file,
      row.names = FALSE
    )

    cat("\n", database, "results for ", label, "\n", sep = "")
    cat(
      "Pathways with adjusted p-value < 0.05:",
      if ("p.adjust" %in% names(result_df)) {
        sum(result_df$p.adjust < 0.05, na.rm = TRUE)
      } else {
        NA
      },
      "\n"
    )

    if (nrow(result_df) > 0) {

      display_cols <- intersect(
        c("ID", "Description", "setSize",
          "NES", "pvalue", "p.adjust"),
        names(result_df)
      )

      result_df <- result_df[
        order(result_df$p.adjust, na.last = TRUE),
        ,
        drop = FALSE
      ]

      print(
        head(result_df[, display_cols, drop = FALSE], 10),
        row.names = FALSE
      )
    }

    invisible(result_df)
  }

  save_and_show(go_bp, "GO_BP")
  save_and_show(kegg, "KEGG")

  # FIXED: Wrapped paths safely in quotes
  saveRDS(
    list(GO_BP = go_bp, KEGG = kegg),
    file.path(output_dir, paste0("GSEA_objects_", label, ".rds"))
  )

  invisible(list(GO_BP = go_bp, KEGG = kegg))
}

# 4. Run the two comparisons
gsea_AD_Old <- run_gsea(
  rank_AD_Old,
  "AD_vs_Old"
)

gsea_Old_Young <- run_gsea(
  rank_Old_Young,
  "Old_vs_Young"
)

# FIXED: Standardized the terminal log message variable name
cat("\nStage 5 completed. Results saved in:", output_dir, "\n")
