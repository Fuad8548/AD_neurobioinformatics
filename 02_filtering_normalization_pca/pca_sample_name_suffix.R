# Stage 2B: Investigate PCA structure by sample-name suffix

library(ggplot2)
library(DESeq2)

pca_df <- plotPCA(
    vsd,
    intgroup = "group",
    returnData = TRUE
)

percent_var <- attr(pca_df, "percentVar")

# Preserve the sample identifiers
pca_df$sample <- rownames(pca_df)

# Extract the trailing A/T code from each sample name
pca_df$sample_code <- sub(
    ".*-\\d+([AT])-(AD|Old|Young)$",
    "\\1",
    pca_df$sample
)

cat("Samples by group and suffix:\n")
print(table(pca_df$group, pca_df$sample_code))

cat("\nSample IDs and suffixes:\n")
print(pca_df[, c("sample", "group", "sample_code")])

# Plot PCA with group colors and A/T shapes
pca_labeled <- ggplot(
    pca_df,
    aes(
        x = PC1,
        y = PC2,
        color = group,
        shape = sample_code,
        label = sample
    )
) +
    geom_point(size = 3.2) +
    geom_text(
        check_overlap = TRUE,
        vjust = -0.8,
        size = 2.6
    ) +
    labs(
        title = "Alzheimer's disease RNA-seq: sample-level PCA",
        subtitle = "Color = group; shape = sample-name suffix",
        x = paste0("PC1: ", percent_var[1], "% variance"),
        y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    theme_bw()

print(pca_labeled)

ggsave(
    file.path(project_dir, "figures", "PCA_sample_labels.png"),
    plot = pca_labeled,
    width = 9,
    height = 6,
    dpi = 300
)
