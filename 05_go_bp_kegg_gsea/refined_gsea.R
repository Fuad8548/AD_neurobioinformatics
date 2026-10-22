# ============================================================
# Stage 6: Refined GSEA with improved p-value resolution
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# FIX: Explicitly define project directory to avoid missing object errors
project_dir <- "results/gsea" 

set.seed(20261021)

gsea_refined <- function(gene_list, label) {

  cat("\nProcessing:", label, "\n")

  go <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    eps = 0,             # Excellent choice for p-value resolution
    by = "fgsea",
    verbose = FALSE
  )

  kegg <- gseKEGG(
    geneList = gene_list,
    organism = "hsa",
    keyType = "ncbi-geneid",
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    eps = 0,             # Excellent choice for p-value resolution
    by = "fgsea",
    verbose = FALSE
  )

  # Print significant pathways from both ends of the ranking
  show_direction <- function(result, database) {

    x <- as.data.frame(result)

    cat("\n", database, ":", label, "\n", sep = "")

    for (direction in c("Positive NES", "Negative NES")) {

      if (direction == "Positive NES") {
        z <- x[!is.na(x$NES) & x$NES > 0 & x$p.adjust < 0.05, , drop = FALSE]
      } else {
        z <- x[!is.na(x$NES) & x$NES < 0 & x$p.adjust < 0.05, , drop = FALSE]
      }

      z <- z[order(z$p.adjust), , drop = FALSE]

      cat("\n", direction, "\n", sep = "")

      if (nrow(z) == 0) {
        cat("No significant pathways in this direction.\n")
      } else {
        cols := intersect(
          c("ID", "Description", "setSize",
            "NES", "pvalue", "p.adjust"),
          names(z)
        )

        print(head(z[, cols, drop = FALSE], 8), row.names = FALSE)
      }
    }
  }

  show_direction(go, "GO_BP")
  show_direction(kegg, "KEGG")

  list(GO_BP = go, KEGG = kegg)
}

# Run the functions (assumes rank_AD_Old and rank_Old_Young exist in memory)
gsea_AD_Old_refined <- gsea_refined(
  rank_AD_Old,
  "AD_vs_Old"
)

gsea_Old_Young_refined <- gsea_refined(
  rank_Old_Young,
  "Old_vs_Young"
)

# Ensure the target directory exists before saving
dir.create(file.path(project_dir, "results"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  list(
    AD_vs_Old = gsea_AD_Old_refined,
    Old_vs_Young = gsea_Old_Young_refined
  ),
  file.path(
    project_dir, "results",
    "GSEA_refined_eps0.rds"
  )
)

cat("\nRefined GSEA finished.\n")
