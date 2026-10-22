# ============================================================
# Alzheimer's Neurobioinformatics
# Stage 2: Filtering, normalization and exploratory PCA
# ============================================================

library(DESeq2)
library(ggplot2)

# 1. Validate the inputs before processing
stopifnot(!anyNA(count_mat))
stopifnot(all(is.finite(count_mat)))
stopifnot(all(count_mat >= 0))
stopifnot(all(count_mat == floor(count_mat)))
stopifnot(!anyDuplicated(gene_symbols))
stopifnot(all(!is.na(gene_symbols) & nzchar(gene_symbols)))

# Assign gene symbols as row names
rownames(count_mat) <- gene_symbols
storage.mode(count_mat) <- "integer"

# 2. Prepare sample metadata
sample_info$group <- factor(
  sample_info$group,
  levels = c("Old", "Young", "AD")
)

coldata <- data.frame(
  group = sample_info$group,
  row.names = sample_info$sample
)

stopifnot(identical(colnames(count_mat), rownames(coldata)))

# 3. Filter low-expression genes
# Retain genes with at least 10 counts in at least 8 samples.
keep <- rowSums(count_mat >= 10) >= 8

cat("Genes before filtering:", nrow(count_mat), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")

count_filtered <- count_mat[keep, , drop = FALSE]

# 4. Create the DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData = coldata,
  design = ~ group
)

# 5. Estimate size factors for library-size normalization
dds <- estimateSizeFactors(dds)

cat("\nDESeq2 size factors:\n")
print(round(sizeFactors(dds), 3))

# Save normalized counts for later inspection
normalized_counts <- counts(dds, normalized = TRUE)

# 6. Variance-stabilizing transformation for exploration
vsd <- vst(dds, blind = TRUE)

# 7. PCA plot
pca_plot <- plotPCA(
  vsd,
  intgroup = "group"
) +
  ggtitle("PCA of Alzheimer's Disease Brain RNA-seq") +
  theme_bw()

print(pca_plot)

# 8. Save the main outputs
dir.create(
  file.path(project_dir, "results"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(project_dir, "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  sample_info,
  file.path(project_dir, "results", "sample_metadata.csv"),
  row.names = FALSE
)

write.csv(
  normalized_counts,
  file.path(project_dir, "results", "normalized_counts.csv")
)

saveRDS(
  dds,
  file.path(project_dir, "results", "dds_preprocessed.rds")
)

saveRDS(
  vsd,
  file.path(project_dir, "results", "vst_object.rds")
)

ggsave(
  file.path(project_dir, "figures", "PCA_initial.png"),
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

cat("\nStage 2 completed. Outputs saved in:", project_dir, "\n")