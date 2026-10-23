# ============================================================
# Stage 9D: APOE4 vs APOE3 cerebral organoid RNA-seq
# Dataset: GSE117588
# ============================================================

library(DESeq2)
library(ggplot2)

if (!exists("project_dir")) project_dir <- getwd()
external_dir <- file.path(project_dir, "data", "external", "GSE117588")

result_dir <- file.path(external_dir, "results")
figure_dir <- file.path(external_dir, "figures")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Prepare the raw count matrix
count_df <- as.data.frame(organoid_raw)
rownames(count_df) <- count_df$geneId
count_df$geneId <- NULL

count_mat <- as.matrix(count_df)
storage.mode(count_mat) <- "numeric"

stopifnot(
    all(is.finite(count_mat)),
    all(count_mat >= 0),
    all(count_mat == round(count_mat))
)

storage.mode(count_mat) <- "integer"

if (anyDuplicated(rownames(count_mat))) {
    stop("Duplicate gene IDs detected. Resolve them before DESeq2.")
}

# 2. Define experimental groups DYNAMICALLY
# Uses the verified metadata from Stage 9B to match the column names perfectly
# ============================================================
# Revised Section 2: Safe Group Assignment without NAs
# ============================================================

if (!exists("sample_metadata")) {
    sample_metadata <- read.csv(file.path(external_dir, "GSE117588_verified_sample_metadata.csv"), stringsAsFactors = FALSE)
}

# Inspect what text is actually inside your metadata column to debug
cat("Raw metadata genotypes found:\n")
print(sample_metadata$apoe_genotype)

# Create a clean group vector based on the metadata content
# If the text contains "E3" or "isogenic", label it E3E3_ISO. Otherwise, E4E4_CTR.
assigned_groups <- ifelse(
    grepl("E3|isogenic|iso", sample_metadata$apoe_genotype, ignore.case = TRUE),
    "E3E3_ISO",
    "E4E4_CTR"
)

sample_info <- data.frame(
    sample = colnames(count_mat),
    group = assigned_groups,
    row.names = colnames(count_mat)
)

# Convert to factor and explicitly set levels (APOE3 is reference)
sample_info$group <- factor(sample_info$group, levels = c("E3E3_ISO", "E4E4_CTR"))

# Safety check: Ensure no NAs slipped through this time
if (any(is.na(sample_info$group))) {
    stop("Critical Error: There are still NA values in your sample groups! Check your metadata file.")
}

stopifnot(
    identical(rownames(sample_info), colnames(count_mat))
)

cat("\nFixed sample groups:\n")
print(sample_info)


# 3. Initial quality control
qc_summary <- data.frame(
    sample = colnames(count_mat),
    group = sample_info$group,
    library_size = colSums(count_mat),
    detected_genes = colSums(count_mat > 0)
)
write.csv(qc_summary, file.path(result_dir, "organoid_QC_summary.csv"), row.names = FALSE)
print(qc_summary)

# 4. Low-expression filtering
keep <- rowSums(count_mat >= 10) >= 3
filtered_counts <- count_mat[keep, , drop = FALSE]

cat("\nGenes before filtering:", nrow(count_mat), "\n")
cat("Genes retained:", nrow(filtered_counts), "\n")

# 5. DESeq2 model
dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = sample_info,
    design = ~group
)
dds <- DESeq(dds)

res_apoe <- results(
    dds,
    contrast = c("group", "E4E4_CTR", "E3E3_ISO"),
    alpha = 0.05
)

res_apoe <- res_apoe[order(res_apoe$padj), ]
res_df <- as.data.frame(res_apoe)
res_df$gene_id <- rownames(res_df)

write.csv(res_df, file.path(result_dir, "DEG_APOE4_vs_APOE3.csv"), row.names = FALSE)

cat("\nDifferential expression summary:\n")
print(summary(res_apoe))

cat("\nSignificant genes at FDR < 0.05:\n")
print(sum(res_apoe$padj < 0.05, na.rm = TRUE))

# 6. Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

saveRDS(dds, file.path(result_dir, "organoid_DESeq2.rds"))
saveRDS(vsd, file.path(result_dir, "organoid_VST.rds"))
saveRDS(res_apoe, file.path(result_dir, "organoid_DEG_results.rds"))

# 7. PCA Plot
pca_df <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = group, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -0.8, show.legend = FALSE) +
    labs(
        title = "APOE genotype and cerebral organoid transcriptomes",
        x = paste0("PC1: ", percent_var[1], "% variance"),
        y = paste0("PC2: ", percent_var[2], "% variance"),
        color = "Genotype"
    ) +
    theme_bw()

print(pca_plot)
ggsave(file.path(figure_dir, "organoid_PCA.png"), pca_plot, width = 7, height = 5, dpi = 300)

# 8. Volcano plot (With fixed scale and explicit coloring setup)
volcano_df <- res_df
volcano_df$significance <- "Not significant"
volcano_df$significance[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & volcano_df$log2FoldChange > 0] <- "Higher in APOE4"
volcano_df$significance[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & volcano_df$log2FoldChange < 0] <- "Higher in APOE3"

volcano_df$neglog10padj <- -log10(pmax(volcano_df$padj, 1e-300))

volcano_plot <- ggplot(volcano_df, aes(log2FoldChange, neglog10padj, color = significance)) +
    geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(0), linetype = 2) + # Added clear center line
    # Added explicit custom colors for clear readability
    scale_color_manual(values = c("Higher in APOE4" = "#d95f02", "Higher in APOE3" = "#7570b3", "Not significant" = "darkgray")) +
    labs(
        title = "APOE4 vs APOE3 organoid differential expression",
        x = "Log2 fold change (APOE4 / APOE3)",
        y = "-Log10 adjusted p-value",
        color = "Expression pattern"
    ) +
    theme_bw()

print(volcano_plot)
ggsave(file.path(figure_dir, "organoid_volcano.png"), volcano_plot, width = 7, height = 5, dpi = 300)

cat("\nStage 9D completed successfully.\n")
