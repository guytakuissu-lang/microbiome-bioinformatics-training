metadata_path <- "data/simulated_ancombc_metadata.tsv"
counts_path <- "data/simulated_ancombc_genus_counts.tsv"
truth_path <- "data/simulated_ancombc_truth.tsv"
report_path <- "results/ancombc_training_taxon_qc.tsv"

metadata <- read.delim(metadata_path, check.names = FALSE)
count_table <- read.delim(counts_path, check.names = FALSE)
truth <- read.delim(truth_path, check.names = FALSE)

required_metadata <- c(
  "sample_id", "group", "age", "sex",
  "metformin", "batch", "library_size"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

if (!"Genus" %in% names(count_table) || !"Genus" %in% names(truth)) {
  stop("Both count and truth tables must contain a Genus column.")
}
if (anyDuplicated(metadata$sample_id)) stop("Duplicated metadata sample IDs.")
if (anyDuplicated(count_table$Genus)) stop("Duplicated count-table genera.")
if (anyDuplicated(truth$Genus)) stop("Duplicated truth-table genera.")

sample_columns <- setdiff(names(count_table), "Genus")
missing_count_samples <- setdiff(metadata$sample_id, sample_columns)
extra_count_samples <- setdiff(sample_columns, metadata$sample_id)
if (length(missing_count_samples) > 0 || length(extra_count_samples) > 0) {
  stop(
    "Sample mismatch. Missing from counts: ",
    paste(missing_count_samples, collapse = ", "),
    "; extra in counts: ", paste(extra_count_samples, collapse = ", ")
  )
}

if (!setequal(count_table$Genus, truth$Genus)) {
  stop("Taxa differ between the count and truth tables.")
}

count_matrix <- as.matrix(count_table[, metadata$sample_id, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- count_table$Genus

if (anyNA(count_matrix)) stop("Missing or non-numeric count values detected.")
if (any(count_matrix < 0)) stop("Negative count values detected.")
if (any(count_matrix %% 1 != 0)) stop("Non-integer count values detected.")

observed_library_sizes <- colSums(count_matrix)
expected_library_sizes <- setNames(metadata$library_size, metadata$sample_id)
if (!all(observed_library_sizes == expected_library_sizes[names(observed_library_sizes)])) {
  stop("Observed counts do not match metadata library sizes.")
}

metadata$group <- factor(metadata$group, levels = c("Control", "T2D"))
metadata$sex <- factor(metadata$sex, levels = c("Female", "Male"))
metadata$metformin <- factor(metadata$metformin, levels = c("No", "Yes"))
metadata$batch <- factor(metadata$batch, levels = c("Batch1", "Batch2"))
if (anyNA(metadata[, c("group", "sex", "metformin", "batch")])) {
  stop("Unexpected categorical level detected in metadata.")
}

design_matrix <- model.matrix(
  ~ group + age + sex + metformin + batch,
  data = metadata
)
design_rank <- qr(design_matrix)$rank
if (design_rank != ncol(design_matrix)) {
  stop("The multivariable design matrix is not full rank.")
}

detection_threshold <- 10
taxon_qc <- data.frame(
  Genus = rownames(count_matrix),
  total_count = rowSums(count_matrix),
  mean_count = rowMeans(count_matrix),
  zero_samples = rowSums(count_matrix == 0),
  samples_with_count_ge_10 = rowSums(count_matrix >= detection_threshold),
  prevalence_ge_10_percent = rowMeans(count_matrix >= detection_threshold) * 100,
  stringsAsFactors = FALSE
)

dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
write.table(
  taxon_qc, report_path, sep = "\t", quote = FALSE, row.names = FALSE
)

cat("======================================\n")
cat("ANCOM-BC2 TRAINING DATA QC\n")
cat("======================================\n\n")
cat("Samples:", nrow(metadata), "\n")
cat("Taxa:", nrow(count_matrix), "\n")
cat("Count cells:", length(count_matrix), "\n")
cat("Zero cells:", sum(count_matrix == 0), "\n")
cat("Library-size range:", min(observed_library_sizes), "to", max(observed_library_sizes), "\n")
cat("Design-matrix rank:", design_rank, "of", ncol(design_matrix), "\n\n")

cat("Group distribution:\n")
print(table(metadata$group))
cat("\nSex by group:\n")
print(table(metadata$group, metadata$sex))
cat("\nBatch by group:\n")
print(table(metadata$group, metadata$batch))
cat("\nMetformin by group:\n")
print(table(metadata$group, metadata$metformin))
cat("\nTaxon QC summary:\n")
print(taxon_qc, row.names = FALSE)
cat("\nAll training-data QC checks passed.\n")
cat("Report saved to:", report_path, "\n")
