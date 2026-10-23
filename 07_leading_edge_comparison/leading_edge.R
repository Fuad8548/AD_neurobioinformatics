# ============================================================
# Stage 7: Leading-edge gene analysis
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Safeguard: Ensure project_dir exists in this session environment
if (!exists("project_dir")) project_dir <- "results/gsea"

# FIX: Align result directory path exactly to match Stage 6 output structures
result_dir <- file.path(project_dir, "results")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# FIX: Read the RDS file generated from the previous session safely
gsea_data_path <- file.path(project_dir, "GSEA_refined_eps0.rds")
if (!file.exists(gsea_data_path)) {
    stop("Missing input file: ", gsea_data_path, "\nDid you run Stage 6 successfully?")
}
gsea_saved_results <- readRDS(gsea_data_path)

# Extract GSEA results out of the loaded RDS structure
go_ad <- gsea_saved_results$AD_vs_Old$GO_BP
go_age <- gsea_saved_results$Old_vs_Young$GO_BP

kegg_ad <- gsea_saved_results$AD_vs_Old$KEGG
kegg_age <- gsea_saved_results$Old_vs_Young$KEGG

# Extract significant pathways occurring in both contrasts with the same direction
shared_same_direction <- function(x, y) {
    dx <- as.data.frame(x)
    dy <- as.data.frame(y)

    if (nrow(dx) == 0 || nrow(dy) == 0) {
        return(character(0))
    }

    ix <- dx$ID[!is.na(dx$p.adjust) & dx$p.adjust < 0.05]
    iy <- dy$ID[!is.na(dy$p.adjust) & dy$p.adjust < 0.05]

    common <- intersect(ix, iy)
    if (length(common) == 0) {
        return(character(0))
    }

    common[vapply(common, function(id) {
        a <- dx[match(id, dx$ID), ]
        b <- dy[match(id, dy$ID), ]
        isTRUE(sign(a$NES) == sign(b$NES))
    }, logical(1))]
}

# Convert an enrichment result's leading-edge IDs to symbols
extract_core <- function(gsea_obj, pathway_id) {
    d <- as.data.frame(gsea_obj)
    if (nrow(d) == 0) {
        return(character(0))
    }

    hit <- d[d$ID == pathway_id, , drop = FALSE]

    if (nrow(hit) == 0 || !("core_enrichment" %in% names(hit)) || is.na(hit$core_enrichment[1])) {
        return(character(0))
    }

    entrez <- unique(strsplit(as.character(hit$core_enrichment[1]), "/", fixed = TRUE)[[1]])

    symbols <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = entrez,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
    )

    unique(unname(symbols[!is.na(symbols)]))
}

# Compare core genes in both contrasts
compare_core <- function(x, y, pathway_ids, database, prefix) {
    dx <- as.data.frame(x)
    ids <- intersect(pathway_ids, shared_same_direction(x, y))

    rows <- list()
    genes <- list()

    if (length(ids) > 0) {
        for (id in ids) {
            a <- extract_core(x, id)
            b <- extract_core(y, id)
            common <- intersect(a, b)
            union_genes <- union(a, b)

            # Fetch description from the data framework safely
            desc <- dx$Description[match(id, dx$ID)]

            rows[[id]] <- data.frame(
                database = database,
                pathway_id = id,
                Description = ifelse(is.na(desc), "Unknown", desc), # FIX: Include Pathway Description
                leading_edge_AD_Old = length(a),
                leading_edge_Old_Young = length(b),
                shared_leading_edge_genes = length(common),
                Jaccard_similarity = if (length(union_genes) > 0) length(common) / length(union_genes) else NA_real_,
                shared_genes = paste(sort(common), collapse = ";"),
                stringsAsFactors = FALSE
            )

            if (length(a) > 0) {
                genes[[paste0(id, "_a")]] <- data.frame(
                    database = database, pathway_id = id, contrast = "AD_vs_Old", SYMBOL = a
                )
            }
            if (length(b) > 0) {
                genes[[paste0(id, "_b")]] <- data.frame(
                    database = database, pathway_id = id, contrast = "Old_vs_Young", SYMBOL = b
                )
            }
        }
    }

    # FIX: Ensure formal column structural generation even if 0 results occur
    summary_df <- if (length(rows) > 0) {
        do.call(rbind, rows)
    } else {
        data.frame(
            database = character(), pathway_id = character(), Description = character(),
            leading_edge_AD_Old = numeric(), leading_edge_Old_Young = numeric(),
            shared_leading_edge_genes = numeric(), Jaccard_similarity = numeric(), shared_genes = character()
        )
    }

    gene_df <- if (length(genes) > 0) {
        do.call(rbind, genes)
    } else {
        data.frame(
            database = character(), pathway_id = character(), contrast = character(), SYMBOL = character()
        )
    }

    write.csv(summary_df, file.path(result_dir, paste0(prefix, "_leading_edge_overlap.csv")), row.names = FALSE)
    write.csv(gene_df, file.path(result_dir, paste0(prefix, "_leading_edge_genes.csv")), row.names = FALSE)

    cat("\n", database, ": leading-edge comparison\n", sep = "")

    if (nrow(summary_df) > 0) {
        print(summary_df[order(-summary_df$shared_leading_edge_genes), ], row.names = FALSE)
    } else {
        cat("No eligible pathways found in both contrasts.\n")
    }

    invisible(summary_df)
}

# 1. Core metabolic pathways
go_core_ids <- c(
    "GO:0006119", # Oxidative phosphorylation
    "GO:0032543", # Mitochondrial translation
    "GO:0042773", # ATP synthesis coupled electron transport
    "GO:0022904" # Respiratory electron transport chain
)

go_core_summary <- compare_core(go_ad, go_age, go_core_ids, "GO_BP", "GO_metabolic")

# 2. KEGG neurodegeneration and energy pathways
kegg_core_ids <- c(
    "hsa00190", # Oxidative phosphorylation
    "hsa05010", # Alzheimer's disease
    "hsa05012", # Parkinson disease
    "hsa05014", # ALS
    "hsa05016", # Huntington disease
    "hsa05020", # Prion disease
    "hsa05022", # Multiple neurodegenerative diseases
    "hsa03050" # Proteasome
)

kegg_core_summary <- compare_core(kegg_ad, kegg_age, kegg_core_ids, "KEGG", "KEGG_neurodegeneration")

cat("\nStage 7 completed.\n")
