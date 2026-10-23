# ============================================================
# Stage 9E: Annotate DEGs with Gene Symbols and Descriptions
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
}

library(org.Hs.eg.db)
library(AnnotationDbi)

# Convert DESeq2 results to a data frame
res_annotated <- as.data.frame(res_apoe)
res_annotated$gene_id <- rownames(res_annotated)

# Safety Check & Fix: Remove Ensembl version extensions (e.g., ENSG00000130203.14 -> ENSG00000130203)
# If version suffixes aren't stripped, database matching will fail completely.
clean_keys <- sub("\\..*$", "", res_annotated$gene_id)

# Optimized Step: Pull all metadata attributes at once using select()
cat("Querying org.Hs.eg.db annotation database...\n")
anno_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = clean_keys,
    columns = c("SYMBOL", "GENENAME", "ENTREZID"),
    keytype = "ENSEMBL"
)

# Resolve multi-mapping items by keeping only the first record for every original Ensembl ID
anno_df <- anno_df[!duplicated(anno_df$ENSEMBL), ]

# Match and merge annotations directly back into the primary dataset
match_idx <- match(clean_keys, anno_df$ENSEMBL)
res_annotated$SYMBOL <- anno_df$SYMBOL[match_idx]
res_annotated$GENENAME <- anno_df$GENENAME[match_idx]
res_annotated$ENTREZID <- anno_df$ENTREZID[match_idx]

# Put the most statistically significant genes first
res_annotated <- res_annotated[order(res_annotated$padj, na.last = TRUE), ]

# Save the annotated results
write.csv(
    res_annotated,
    file.path(external_dir, "results", "DEG_APOE4_vs_APOE3_annotated.csv"),
    row.names = FALSE
)

# Show the top 25 results in R console
print(
    head(
        res_annotated[, c("gene_id", "SYMBOL", "GENENAME", "baseMean", "log2FoldChange", "padj")],
        25
    )
)

# Check annotation coverage metrics
mapped_count <- sum(!is.na(res_annotated$SYMBOL))
cat(
    "\nAnnotation Verification Metrics:\n",
    "Genes with successfully mapped symbols:", mapped_count, "out of", nrow(res_annotated),
    paste0("(", round((mapped_count / nrow(res_annotated)) * 100, 1), "% coverage)\n")
)
