# Install packages if needed
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
    BiocManager::install("enrichplot")
}

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# Create a ranked list using DESeq2's Wald statistic
gsea_df <- as.data.frame(res_apoe)
gsea_df$ENSEMBL <- rownames(gsea_df)

# Map Ensembl IDs to Entrez IDs
gsea_df$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = gsea_df$ENSEMBL,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

gsea_df <- gsea_df[
    !is.na(gsea_df$ENTREZID) &
        !is.na(gsea_df$stat),
]

# Keep one entry per Entrez ID
gsea_df <- gsea_df[
    order(abs(gsea_df$stat), decreasing = TRUE),
]
gsea_df <- gsea_df[
    !duplicated(gsea_df$ENTREZID),
]

gene_rank <- gsea_df$stat
names(gene_rank) <- gsea_df$ENTREZID
gene_rank <- sort(gene_rank, decreasing = TRUE)

cat("Genes in ranked list:", length(gene_rank), "\n")

# Gene Ontology Biological Process GSEA
gsea_go <- gseGO(
    geneList = gene_rank,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
)

# KEGG pathway GSEA
gsea_kegg <- gseKEGG(
    geneList = gene_rank,
    organism = "hsa",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
)

# Create results folder if it does not exist
dir.create(
    file.path(external_dir, "results"),
    recursive = TRUE,
    showWarnings = FALSE
)
dir.create(
    file.path(external_dir, "figures"),
    recursive = TRUE,
    showWarnings = FALSE
)

# Save pathway tables
write.csv(
    as.data.frame(gsea_go),
    file.path(
        external_dir, "results",
        "GSEA_GO_Biological_Process.csv"
    ),
    row.names = FALSE
)

write.csv(
    as.data.frame(gsea_kegg),
    file.path(
        external_dir, "results",
        "GSEA_KEGG.csv"
    ),
    row.names = FALSE
)

# Display the top enriched pathways
go_table <- as.data.frame(gsea_go)
kegg_table <- as.data.frame(gsea_kegg)

print(head(
    go_table[, intersect(
        c("ID", "Description", "NES", "p.adjust", "core_enrichment"),
        names(go_table)
    )],
    15
))

print(head(
    kegg_table[, intersect(
        c("ID", "Description", "NES", "p.adjust", "core_enrichment"),
        names(kegg_table)
    )],
    15
))

# Plot the top GO pathways
if (nrow(go_table) > 0) {
    p_go <- dotplot(gsea_go, showCategory = 15) +
        ggtitle("GSEA: GO Biological Processes")

    ggsave(
        file.path(external_dir, "figures", "GSEA_GO_BP.png"),
        p_go,
        width = 10, height = 7, dpi = 300
    )
}

# Plot the top KEGG pathways
if (nrow(kegg_table) > 0) {
    p_kegg <- dotplot(gsea_kegg, showCategory = 15) +
        ggtitle("GSEA: KEGG Pathways")

    ggsave(
        file.path(external_dir, "figures", "GSEA_KEGG.png"),
        p_kegg,
        width = 10, height = 7, dpi = 300
    )
}
