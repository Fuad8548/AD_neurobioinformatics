# ============================================================
# Stage 9F: Leading-Edge Gene Validation & Pathway Visualization
# ============================================================

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(DESeq2) # Required for reading the dds object and counts() function

# 1. Re-establish environmental directories safely
if (!exists("project_dir")) project_dir <- getwd()
external_dir <- file.path(project_dir, "data", "external", "GSE117588")
result_dir <- file.path(external_dir, "results")
figure_dir <- file.path(external_dir, "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Paths to stored RDS objects
gsea_go_rds <- file.path(result_dir, "GSEA_GO_BP_object.rds")
annotated_deg_rds <- file.path(result_dir, "organoid_DEG_annotated.rds")
deseq2_dds_rds <- file.path(result_dir, "organoid_DESeq2.rds")

# Verification and loading step
if (!all(file.exists(c(gsea_go_rds, annotated_deg_rds, deseq2_dds_rds)))) {
    stop("Missing required RDS files in results directory. Please ensure previous stages ran completely.")
}

cat("Loading saved analysis objects from RDS files...\n")
gsea_go <- readRDS(gsea_go_rds)
res_annotated <- readRDS(annotated_deg_rds)
dds <- readRDS(deseq2_dds_rds)

# ============================================================
# Part A: Extract Leading-Edge Genes from Core Targets
# ============================================================

# Configuration updated to match your exact structural pathway ID (GO:0043062)
selected_go <- c(
    "GO:0097529", # Myeloid leukocyte migration
    "GO:0043062", # Extracellular structure organization (Fixed ID)
    "GO:0099504", # Synaptic vesicle cycle
    "GO:0007416" # Synapse assembly
)

go_df <- as.data.frame(gsea_go)
selected_results <- go_df[
    go_df$ID %in% selected_go,
    c("ID", "Description", "NES", "p.adjust", "core_enrichment")
]

cat("\nExtracting leading-edge metadata for selected categories:\n")
print(selected_results[, c("ID", "Description", "NES")])

# Loop to map Entrez IDs back to clean symbols
leading_edge_list <- lapply(
    seq_len(nrow(selected_results)),
    function(i) {
        if (is.na(selected_results$core_enrichment[i]) || selected_results$core_enrichment[i] == "") {
            return(NULL)
        }

        ids <- unlist(strsplit(selected_results$core_enrichment[i], "/"))
        symbols <- mapIds(
            org.Hs.eg.db,
            keys = ids,
            column = "SYMBOL",
            keytype = "ENTREZID",
            multiVals = "first"
        )

        data.frame(
            pathway = selected_results$Description[i],
            ENTREZID = ids,
            SYMBOL = unname(symbols),
            stringsAsFactors = FALSE
        )
    }
)

leading_edge_list <- leading_edge_list[!sapply(leading_edge_list, is.null)]
leading_edge_table <- do.call(rbind, leading_edge_list)

write.csv(
    leading_edge_table,
    file.path(result_dir, "selected_pathway_leading_edge_genes.csv"),
    row.names = FALSE
)

# ============================================================
# Part B: Print and Validate Normalized Expression Counts
# ============================================================

genes_to_check <- c("H19", "NLRP2", "REST")
gene_ids <- res_annotated$gene_id[match(genes_to_check, res_annotated$SYMBOL)]
gene_ids <- gene_ids[!is.na(gene_ids)] # Remove unmapped IDs

# Extract normalized matrix counts directly from the restored dds model
dds_counts <- counts(dds, normalized = TRUE)
valid_ids <- gene_ids[gene_ids %in% rownames(dds_counts)]

if (length(valid_ids) > 0) {
    count_check <- dds_counts[valid_ids, , drop = FALSE]

    # Change row names back to human gene symbols for cleaner terminal display
    rownames(count_check) <- res_annotated$SYMBOL[match(rownames(count_check), res_annotated$gene_id)]

    cat("\nValidated Target Gene Fingerprints:\n")
    print(round(count_check, 1))
} else {
    warning("None of the specified genes could be identified in the expression matrix scale.")
}

# ============================================================
# Part C: Generate and Save Publication-Grade GSEA Mountain Plots
# ============================================================

cat("\nGenerating running-sum enrichment graphs...\n")
available_ids <- go_df$ID

pathway_plots_config <- data.frame(
    id = c("GO:0099504", "GO:0097529", "GO:0043062"),
    title = c("Synaptic vesicle cycle", "Myeloid leukocyte migration", "Extracellular structure organization"),
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(pathway_plots_config))) {
    pathway_id <- pathway_plots_config$id[i]
    pathway_title <- pathway_plots_config$title[i]

    if (pathway_id %in% available_ids) {
        p <- gseaplot2(
            gsea_go,
            geneSetID    = pathway_id,
            title        = pathway_title,
            pvalue_table = TRUE
        )

        safe_name <- gsub("[^A-Za-z0-9]+", "_", pathway_title)

        ggsave(
            filename = file.path(figure_dir, paste0("GSEA_", safe_name, ".png")),
            plot     = p,
            width    = 9,
            height   = 6,
            dpi      = 300
        )
    }
}

cat("\nAll tasks completed successfully. Selected matrices and plots saved.\n")
