
# Stage 1B: Validate sample metadata and count matrix

# 1. Extract sample names
sample_names <- names(counts_raw)[-1]

# 2. Infer groups from the sample-name suffix
group <- sub(".*-(AD|Old|Young)$", "\\1", sample_names)

# Create a metadata table
sample_info <- data.frame(
  sample = sample_names,
  group = group,
  stringsAsFactors = FALSE
)

cat("Sample groups:\n")
print(table(sample_info$group))

cat("\nAll sample assignments:\n")
print(sample_info, row.names = FALSE)

# 3. Convert the count columns to a numeric matrix
count_df <- counts_raw[, -1, drop = FALSE]

count_mat <- as.matrix(
  data.frame(
    lapply(count_df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
)

colnames(count_mat) <- names(count_df)

gene_symbols <- trimws(as.character(counts_raw$refGene))

cat("\nCount matrix dimensions:\n")
print(dim(count_mat))

cat("\nMissing count values:", sum(is.na(count_mat)), "\n")

cat(
  "Non-integer count values:",
  sum(count_mat[!is.na(count_mat)] !=
        floor(count_mat[!is.na(count_mat)])),
  "\n"
)

cat(
  "Negative count values:",
  sum(count_mat < 0, na.rm = TRUE),
  "\n"
)

# 4. Check gene identifiers
cat(
  "\nMissing or blank gene symbols:",
  sum(is.na(gene_symbols) | gene_symbols == ""),
  "\n"
)

valid_symbols <- gene_symbols[
  !is.na(gene_symbols) & gene_symbols != ""
]

symbol_frequency <- sort(
  table(valid_symbols),
  decreasing = TRUE
)

cat(
  "Repeated gene-symbol rows:",
  sum(duplicated(valid_symbols)),
  "\n"
)

cat("\nGene symbols occurring more than once:\n")
print(head(symbol_frequency[symbol_frequency > 1], 15))

# 5. Inspect library sizes
library_sizes <- colSums(count_mat, na.rm = TRUE)

cat("\nLibrary sizes:\n")
print(sort(library_sizes))

cat("\nSummary of library sizes:\n")
print(summary(library_sizes))


