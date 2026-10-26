# ============================================================
# Stage 9B: Verify APOE genotype for every organoid sample
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

library(GEOquery)

# 1. Fast download using GSEMatrix = TRUE
cat("Downloading pre-parsed GEO metadata...\n")
gse_list <- getGEO("GSE117588", GSEMatrix = TRUE)
gse_meta <- gse_list[[1]]
all_meta <- pData(gse_meta)

# 2. Extract and format target fields dynamically
characteristic_cols <- grep("^characteristics_ch1", names(all_meta), value = TRUE)

sample_metadata <- data.frame(
    accession = as.character(all_meta$geo_accession),
    sample = as.character(all_meta$title),
    source_name = as.character(all_meta$source_name_ch1),
    characteristics = apply(all_meta[characteristic_cols], 1, function(x) {
        paste(na.omit(x), collapse = "; ")
    }),
    stringsAsFactors = FALSE
)

# 3. Match against the count matrix columns from Stage 9A
# Ensure organoid_raw exists from your previous script execution
if (!exists("organoid_raw")) {
    stop("Error: 'organoid_raw' matrix not found. Please run Stage 9A first.")
}

count_samples <- names(organoid_raw)[-1]

# Safety Check: Verify if columns match titles or accessions
if (!any(count_samples %in% sample_metadata$sample) && any(count_samples %in% sample_metadata$accession)) {
    cat("Notice: Count matrix uses GSM accession IDs instead of sample titles. Adjusting match...\n")
    sample_metadata$match_key <- sample_metadata$accession
} else {
    sample_metadata$match_key <- sample_metadata$sample
}

# Align metadata perfectly with the count matrix order
match_indices <- match(count_samples, sample_metadata$match_key)

# Warning if any sample is missing metadata
if (any(is.na(match_indices))) {
    warning("Some samples in the count matrix do not have matching metadata in GEO!")
}

sample_metadata <- sample_metadata[match_indices, c("accession", "sample", "source_name", "characteristics"), drop = FALSE]

# 4. Save and inspect
cat("\nVerified sample annotations:\n")
print(head(sample_metadata), row.names = FALSE)

write.csv(
    sample_metadata,
    file.path(external_dir, "GSE117588_verified_sample_metadata.csv"),
    row.names = FALSE
)
cat("\nMetadata successfully saved.\n")
