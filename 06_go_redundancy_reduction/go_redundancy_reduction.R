# ============================================================
# Stage 6: GO redundancy reduction and pathway convergence
# ============================================================

library(clusterProfiler)
library(GOSemSim) # Required for the backend semantic similarity metrics

# Safeguard: Ensure project_dir exists in this session environment
if (!exists("project_dir")) project_dir <- "results/gsea"

# FIX: Align result directory logic with where the data was actually saved
result_dir <- file.path(project_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# FIX: Read the RDS file generated from the previous script
gsea_data_path <- file.path(project_dir, "GSEA_refined_eps0.rds")
if (!file.exists(gsea_data_path)) {
    stop("Missing input file: ", gsea_data_path, "\nDid you run the previous stage successfully?")
}
gsea_saved_results <- readRDS(gsea_data_path)

# 1. Retrieve the existing GSEA results safely from the loaded file
go_ad <- gsea_saved_results$AD_vs_Old$GO_BP
go_age <- gsea_saved_results$Old_vs_Young$GO_BP

kegg_ad <- gsea_saved_results$AD_vs_Old$KEGG
kegg_age <- gsea_saved_results$Old_vs_Young$KEGG

# 2. Reduce redundant GO Biological Process terms safely
simplify_gsea <- function(gsea_obj) {
    tryCatch(
        simplify(
            gsea_obj,
            cutoff = 0.7,
            by = "p.adjust",
            select_fun = min,
            measure = "Wang"
        ),
        error = function(e) {
            message("Simplification skipped or failed for this object: ", conditionMessage(e))
            return(gsea_obj)
        }
    )
}

go_ad_simple <- simplify_gsea(go_ad)
go_age_simple <- simplify_gsea(go_age)

write.csv(
    as.data.frame(go_ad_simple),
    file.path(result_dir, "GO_BP_AD_vs_Old_simplified.csv"),
    row.names = FALSE
)

write.csv(
    as.data.frame(go_age_simple),
    file.path(result_dir, "GO_BP_Old_vs_Young_simplified.csv"),
    row.names = FALSE
)

# 3. Display the strongest significant terms after simplification
show_top <- function(x, label, n = 10) {
    d <- as.data.frame(x)
    if (nrow(d) == 0) {
        return(cat("\n", label, ": No pathways found.\n", sep = ""))
    }

    d <- d[!is.na(d$p.adjust) & d$p.adjust < 0.05, ]
    if (nrow(d) == 0) {
        return(cat("\n", label, ": No significant pathways found.\n", sep = ""))
    }

    d <- d[order(d$p.adjust), ]

    cat("\n", label, ": significant nonredundant GO terms\n", sep = "")
    cols_to_show <- intersect(c("ID", "Description", "NES", "p.adjust"), names(d))
    print(head(d[, cols_to_show], n), row.names = FALSE)
}

show_top(go_ad_simple, "AD vs older controls")
show_top(go_age_simple, "Older vs younger controls")

# 4. Prepare significant pathway tables
sig_table <- function(x, suffix) {
    d <- as.data.frame(x)
    if (nrow(d) == 0) {
        return(data.frame(ID = character(), Description = character(), NES = numeric(), p.adjust = numeric()))
    }

    d <- d[
        !is.na(d$p.adjust) & d$p.adjust < 0.05,
        c("ID", "Description", "NES", "p.adjust"),
        drop = FALSE
    ]

    if (nrow(d) == 0) {
        return(data.frame(ID = character(), Description = character(), NES = numeric(), p.adjust = numeric()))
    }

    names(d)[names(d) %in% c("NES", "p.adjust")] <- paste0(names(d)[names(d) %in% c("NES", "p.adjust")], suffix)
    d
}

# Compare shared GO terms using both ID and Description keys
go_ad_sig <- sig_table(go_ad, "_AD_Old")
go_age_sig <- sig_table(go_age, "_Old_Young")

go_overlap <- merge(go_ad_sig, go_age_sig, by = c("ID", "Description"))

# FIX: Safeguard empty row operations before adding columns or sorting
if (nrow(go_overlap) > 0) {
    go_overlap$direction_pattern <- ifelse(
        sign(go_overlap$NES_AD_Old) == sign(go_overlap$NES_Old_Young),
        "Same NES direction",
        "Opposite NES directions"
    )
    go_overlap <- go_overlap[order(go_overlap$p.adjust_AD_Old), ]
}

write.csv(
    go_overlap,
    file.path(result_dir, "GO_pathway_convergence.csv"),
    row.names = FALSE
)

# 5. Repeat for shared KEGG pathways
kegg_ad_sig <- sig_table(kegg_ad, "_AD_Old")
kegg_age_sig <- sig_table(kegg_age, "_Old_Young")

kegg_overlap <- merge(kegg_ad_sig, kegg_age_sig, by = c("ID", "Description"))

# FIX: Safeguard empty row operations
if (nrow(kegg_overlap) > 0) {
    kegg_overlap$direction_pattern <- ifelse(
        sign(kegg_overlap$NES_AD_Old) == sign(kegg_overlap$NES_Old_Young),
        "Same NES direction",
        "Opposite NES directions"
    )
    kegg_overlap <- kegg_overlap[order(kegg_overlap$p.adjust_AD_Old), ]
}

write.csv(
    kegg_overlap,
    file.path(result_dir, "KEGG_pathway_convergence.csv"),
    row.names = FALSE
)

cat("\nShared significant GO pathways:", nrow(go_overlap), "\n")
cat("Shared significant KEGG pathways:", nrow(kegg_overlap), "\n")

if (nrow(go_overlap) > 0) {
    cat("\nTop shared GO pathways:\n")
    print(head(go_overlap, 10), row.names = FALSE)
}

if (nrow(kegg_overlap) > 0) {
    cat("\nTop shared KEGG pathways:\n")
    print(head(kegg_overlap, 10), row.names = FALSE)
}

cat("\nStage 6 completed. Files saved in:", result_dir, "\n")
