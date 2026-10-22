# Stage 1: Download and inspect the raw gene-count matrix

project_dir <- "AD_Neurobioinformatics"

dir.create(
  file.path(project_dir, "data", "raw"),
  recursive = TRUE,
  showWarnings = FALSE
)

count_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE159nnn/GSE159699/suppl/",
  "GSE159699_summary_count.star.txt.gz"
)

count_file <- file.path(
  project_dir, "data", "raw",
  "GSE159699_summary_count.star.txt.gz"
)

if (!file.exists(count_file)) {
  download.file(
    count_url,
    destfile = count_file,
    mode = "wb"
  )
}

# Read the compressed count matrix
con <- gzfile(count_file, open = "rt")

counts_raw <- read.delim(
  con,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

close(con)

# Inspect the structure
cat("Matrix dimensions:", dim(counts_raw), "\n\n")

cat("First 10 column names:\n")
print(head(names(counts_raw), 10))

cat("\nLast 5 column names:\n")
print(tail(names(counts_raw), 5))

cat("\nFirst 5 rows and columns:\n")
print(head(counts_raw[, seq_len(min(6, ncol(counts_raw)))]))


