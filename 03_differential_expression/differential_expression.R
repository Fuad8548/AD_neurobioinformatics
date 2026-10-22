# ============================================================
# Project 2: Alzheimer's Neurobioinformatics
# Stage 3: Primary differential expression
# Contrast: Alzheimer's disease vs older controls
# ============================================================

library(DESeq2)

dds <- readRDS("results/dds_preprocessed.rds")

# 1. Fit the differential-expression model
dds <- DESeq(dds)

# 2. Extract the AD vs Old contrast
res_AD_Old <- results(
    dds,
    contrast = c("group", "AD", "Old"),
    alpha = 0.05
)

# 3. Summarize the results
cat("\nDifferential-expression summary:\n")
print(summary(res_AD_Old))

cat(
    "\nGenes with adjusted p-value < 0.05:",
    sum(res_AD_Old$padj < 0.05, na.rm = TRUE),
    "\n"
)

# 4. Organize and rank the results
res_AD_Old_df <- data.frame(
    gene = rownames(res_AD_Old),
    as.data.frame(res_AD_Old),
    check.names = FALSE
)

res_AD_Old_df <- res_AD_Old_df[
    order(res_AD_Old_df$padj, na.last = TRUE),
]

# 5. Inspect the top results
cat("\nTop 15 genes by adjusted p-value:\n")

print(
    head(
        res_AD_Old_df[
            , c(
                "gene", "baseMean", "log2FoldChange",
                "lfcSE", "stat", "pvalue", "padj"
            )
        ],
        15
    )
)

# 6. Save the results
write.csv(
    res_AD_Old_df,
    "results/DESeq2_AD_vs_Old.csv",
    row.names = FALSE
)

saveRDS(
    dds,
    "results/dds_fitted.rds"
)

# 7. Save an MA plot
png(
    "figures/MA_AD_vs_Old.png",
    width = 1800,
    height = 1400,
    res = 220
)

plotMA(
    res_AD_Old,
    ylim = c(-5, 5),
    main = "Alzheimer's Disease vs Older Controls"
)

dev.off()

cat("\nStage 3 completed.\n")
