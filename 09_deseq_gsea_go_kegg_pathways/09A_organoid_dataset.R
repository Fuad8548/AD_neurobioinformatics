# ============================================================
# Stage 9A: Inspect the professor-aligned organoid dataset
# GEO: GSE117588
# ============================================================

# Ensure project_dir is defined (e.g., project_dir <- ".")
if (!exists("project_dir")) project_dir <- "data/external"

external_dir <- file.path(project_dir, "data", "external", "GSE117588")
dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)

organoid_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    "GSE117nnn/GSE117588/suppl/",
    "GSE117588_organoids-gene-counts.tsv.gz"
)

organoid_file <- file.path(external_dir, "GSE117588_organoids-gene-counts.tsv.gz")

if (!file.exists(organoid_file)) {
    download.file(organoid_url, destfile = organoid_file, mode = "wb")
}

# R automatically handles .gz decompression directly from the file path
organoid_raw <- read.delim(
    organoid_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("Matrix dimensions:\n")
print(dim(organoid_raw))

cat("\nColumn names (first 10):\n")
print(head(names(organoid_raw), 10))

cat("\nFirst six rows:\n")
print(head(organoid_raw[, 1:min(5, ncol(organoid_raw))])) # Safe printing for large tables

cat("\nData types summary:\n")
print(table(sapply(organoid_raw, class))) # Summarizes column types instead of listing all

cat("\nInspection completed.\n")
