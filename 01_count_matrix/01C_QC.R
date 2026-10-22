
# Stage 1C: Final basic QC checks

cat("Sample counts by group:\n")
print(table(sample_info$group))

cat("\nSample order matches count matrix: ",
    identical(colnames(count_mat), sample_info$sample), "\n")

cat("Missing count values: ",
    sum(is.na(count_mat)), "\n")

cat("Non-finite count values: ",
    sum(!is.finite(count_mat)), "\n")

cat("Negative count values: ",
    sum(count_mat < 0, na.rm = TRUE), "\n")

cat("Non-integer count values: ",
    sum(count_mat[is.finite(count_mat)] !=
        floor(count_mat[is.finite(count_mat)])), "\n")

cat("Missing or blank gene symbols: ",
    sum(is.na(gene_symbols) | trimws(gene_symbols) == ""), "\n")

# Inspect duplicate gene symbols
valid_symbols <- gene_symbols[
  !is.na(gene_symbols) & trimws(gene_symbols) != ""
]

symbol_counts <- sort(table(valid_symbols), decreasing = TRUE)

cat("Number of repeated gene-symbol rows: ",
    sum(duplicated(valid_symbols)), "\n")

cat("\nExamples of repeated symbols:\n")
print(head(symbol_counts[symbol_counts > 1], 15))




















