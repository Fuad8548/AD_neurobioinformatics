# ============================================================
# Stage 8: Integrate mitochondrial leading-edge gene results
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# Safeguard: Ensure project_dir exists in this session environment
if (!exists("project_dir")) project_dir <- "results/gsea"

# FIX: Set up paths and ensure BOTH output folders are created to avoid ggsave crashes
result_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# FIX: Load GSEA results from Stage 6/7 RDS file
gsea_data_path <- file.path(project_dir, "GSEA_refined_eps0.rds")
if (!file.exists(gsea_data_path)) {
    stop("Missing input file: ", gsea_data_path, "\nDid you run previous stages successfully?")
}
gsea_saved_results <- readRDS(gsea_data_path)

# Extract objects out of the loaded GSEA structure safely
go_ad <- gsea_saved_results$AD_vs_Old$GO_BP
gsea_AD_Old_refined <- gsea_saved_results$AD_vs_Old
gsea_Old_Young_refined <- gsea_saved_results$Old_vs_Young

# FIX: Load your raw differential expression tables from memory/disk
# Update these paths to target your actual DESeq2/EdgeR result filenames
if (!file.exists("res_AD_Old.rds") || !file.exists("res_Old_Young.rds")) {
    stop("Missing raw Differential Expression RDS files ('res_AD_Old.rds' / 'res_Old_Young.rds') required for gene expression joins!")
}
res_AD_Old <- readRDS("res_AD_Old.rds")
res_Old_Young <- readRDS("res_Old_Young.rds")

# 1. Extract the leading-edge genes from GO respiratory electron transport
go_results <- as.data.frame(go_ad)
pathway_id <- "GO:0022904"
pathway_row <- go_results[go_results$ID == pathway_id, , drop = FALSE]

if (nrow(pathway_row) != 1) {
    stop("Pathway ID ", pathway_id, " was not found or is not significant in the GO object.")
}

core_entrez <- unique(strsplit(pathway_row$core_enrichment[1], "/", fixed = TRUE)[[1]])

core_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = core_entrez,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
)
core_symbols <- unique(unname(core_symbols[!is.na(core_symbols)]))

# 2. Prepare gene-level results for each comparison
make_de_table <- function(res, suffix) {
    x <- as.data.frame(res)

    # Ensure column structural robustness depending on engine tool types
    df <- data.frame(
        SYMBOL = rownames(x),
        log2FC = x$log2FoldChange,
        Wald_stat = if ("stat" %in% names(x)) x$stat else NA_real_,
        padj = x$padj,
        stringsAsFactors = FALSE
    )

    names(df) <- c("SYMBOL", paste0("log2FC_", suffix), paste0("Wald_stat_", suffix), paste0("padj_", suffix))
    return(df)
}

ad_table <- make_de_table(res_AD_Old, "AD_Old")
age_table <- make_de_table(res_Old_Young, "Old_Young")

oxphos_table <- merge(data.frame(SYMBOL = core_symbols), ad_table, by = "SYMBOL", all.x = TRUE)
oxphos_table <- merge(oxphos_table, age_table, by = "SYMBOL", all.x = TRUE)

# 3. Classify the observed patterns safely accounting for NA cases
sig_ad <- !is.na(oxphos_table$padj_AD_Old) & oxphos_table$padj_AD_Old < 0.05
sig_age <- !is.na(oxphos_table$padj_Old_Young) & oxphos_table$padj_Old_Young < 0.05

# FIX: Safeguard sign evaluation from unexpected NAs crashing ifelse logic
same_direction <- !is.na(oxphos_table$log2FC_AD_Old) & !is.na(oxphos_table$log2FC_Old_Young) &
    (sign(oxphos_table$log2FC_AD_Old) == sign(oxphos_table$log2FC_Old_Young))

oxphos_table$pattern <- ifelse(
    sig_ad & sig_age & same_direction, "Both comparisons, same direction",
    ifelse(sig_ad & sig_age & !same_direction, "Both comparisons, opposite directions",
        ifelse(sig_ad & !sig_age, "AD comparison only",
            ifelse(!sig_ad & sig_age, "Aging comparison only", "Neither significant at FDR 0.05")
        )
    )
)

# 4. Save the integrated results
oxphos_table <- oxphos_table[order(oxphos_table$padj_AD_Old, na.last = TRUE), ]
write.csv(oxphos_table, file.path(result_dir, "OXPHOS_gene_level_integration.csv"), row.names = FALSE)

cat("Number of leading-edge genes:", nrow(oxphos_table), "\n\n")
cat("Patterns across both comparisons:\n")
print(table(oxphos_table$pattern))

cat("\nTop 15 genes by AD-vs-older-control adjusted p-value:\n")
print(head(oxphos_table, 15), row.names = FALSE)

# 5. Compare effect sizes across contrasts
valid <- complete.cases(oxphos_table$log2FC_AD_Old, oxphos_table$log2FC_Old_Young)

if (sum(valid) >= 3) {
    cat("\nSpearman correlation between the two contrasts:\n")
    print(cor.test(oxphos_table$log2FC_AD_Old[valid], oxphos_table$log2FC_Old_Young[valid], method = "spearman", exact = FALSE))
}

effect_plot <- ggplot(oxphos_table, aes(x = log2FC_Old_Young, y = log2FC_AD_Old, color = pattern)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(size = 2.5, alpha = 0.85) +
    labs(
        title = "Shared mitochondrial genes across disease and aging",
        x = "Log2 fold change: older vs younger controls",
        y = "Log2 fold change: Alzheimer's disease vs older controls",
        color = "Observed pattern"
    ) +
    theme_bw()

print(effect_plot)
ggsave(file.path(figure_dir, "OXPHOS_effect_comparison.png"), plot = effect_plot, width = 8, height = 6, dpi = 300)

# 6. Check the REST gene's own expression, if available safely
res_ad_df <- as.data.frame(res_AD_Old)
res_age_df <- as.data.frame(res_Old_Young)

if ("REST" %in% rownames(res_ad_df) && "REST" %in% rownames(res_age_df)) {
    rest_table <- rbind(
        data.frame(contrast = "AD vs older controls", SYMBOL = "REST", res_ad_df["REST", , drop = FALSE]),
        data.frame(contrast = "Older vs younger controls", SYMBOL = "REST", res_age_df["REST", , drop = FALSE])
    )
    cat("\nREST gene-level differential expression:\n")
    print(rest_table, row.names = FALSE)
} else {
    cat("\nREST was not found among the count-matrix row names.\n")
}

# 7. Inspect neural pathways
show_neural_terms <- function(gsea_obj, label) {
    d <- as.data.frame(gsea_obj)
    if (nrow(d) == 0) {
        cat("\nNeural-process pathways:", label, "\nNo pathways found.\n")
        return(invisible(NULL))
    }

    keep <- grepl("neuron|synap|neurogen|axon|dendrit|trans-synaptic", d$Description, ignore.case = TRUE)
    d <- d[keep, , drop = FALSE]

    cat("\nNeural-process pathways:", label, "\n")
    if (nrow(d) == 0) {
        cat("No matching pathways were returned.\n")
        return(invisible(NULL))
    }

    d <- d[order(d$p.adjust, na.last = TRUE), ]
    cols <- intersect(c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust"), names(d))
    print(head(d[, cols, drop = FALSE], 20), row.names = FALSE)
    invisible(d)
}

show_neural_terms(gsea_AD_Old_refined$GO_BP, "AD vs older controls")
show_neural_terms(gsea_Old_Young_refined$GO_BP, "Older vs younger controls")

cat("\nStage 8 completed.\n")
