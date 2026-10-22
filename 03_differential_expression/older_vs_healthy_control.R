# ============================================================
# Stage 4: Normal-aging comparison
# Older healthy controls vs younger healthy controls
# ============================================================

library(DESeq2)

# 1. Run the secondary contrast
res_Old_Young <- results(
  dds,
  contrast = c("group", "Old", "Young"),
  alpha = 0.05
)

cat("\nOlder controls vs younger controls:\n")
print(summary(res_Old_Young))

cat(
  "\nSignificant genes (adjusted p < 0.05):",
  sum(res_Old_Young$padj < 0.05, na.rm = TRUE),
  "\n"
)

# 2. Organize the results
res_Old_Young_df <- as.data.frame(res_Old_Young)
res_Old_Young_df$gene <- rownames(res_Old_Young_df)

res_Old_Young_df <- res_Old_Young_df[
  order(
    res_Old_Young_df$padj,
    -abs(res_Old_Young_df$stat),
    na.last = TRUE
  ),
]

# 3. Inspect the top 15 genes
cat("\nTop 15 genes: older vs younger controls\n")

print(
  head(
    res_Old_Young_df[
      , c(
        "gene", "baseMean", "log2FoldChange",
        "lfcSE", "stat", "pvalue", "padj"
      )
    ],
    15
  ),
  row.names = FALSE
)

# 4. Save the complete results
write.csv(
  res_Old_Young_df,
  file.path(
    project_dir,
    "results",
    "DESeq2_Old_vs_Young.csv"
  ),
  row.names = FALSE
)

# 5. Save a consistent MA plot
png(
  file.path(
    project_dir,
    "figures",
    "MA_Old_vs_Young.png"
  ),
  width = 1800,
  height = 1400,
  res = 220
)

plotMA(
  res_Old_Young,
  alpha = 0.05,
  ylim = c(-5, 5),
  main = "Older vs Younger Healthy Controls"
)

dev.off()

cat("\nStage 4 completed.\n")