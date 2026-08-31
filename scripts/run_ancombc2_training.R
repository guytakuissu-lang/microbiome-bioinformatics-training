suppressPackageStartupMessages(library(ANCOMBC))

counts_path <- "data/simulated_ancombc_genus_counts.tsv"
metadata_path <- "data/simulated_ancombc_metadata.tsv"
truth_path <- "data/simulated_ancombc_truth.tsv"

filter_path <- "results/ancombc2_training_filter.tsv"
full_results_path <- "results/ancombc2_training_full_results.tsv"
group_results_path <- "results/ancombc2_training_group_results.tsv"

detection_threshold <- 10
minimum_detection_prevalence <- 0.20

count_table <- read.delim(counts_path, check.names = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE)
truth <- read.delim(truth_path, check.names = FALSE)

if (!"Genus" %in% names(count_table)) {
  stop("The count table does not contain a Genus column.")
}

count_matrix <- as.matrix(count_table[, -1, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- count_table$Genus

if (!identical(colnames(count_matrix), metadata$sample_id)) {
  stop("Count-table columns and metadata sample IDs are not identically ordered.")
}

rownames(metadata) <- metadata$sample_id
metadata$group <- factor(metadata$group, levels = c("Control", "T2D"))
metadata$sex <- factor(metadata$sex, levels = c("Female", "Male"))
metadata$metformin <- factor(metadata$metformin, levels = c("No", "Yes"))
metadata$batch <- factor(metadata$batch, levels = c("Batch1", "Batch2"))
metadata$age_centered <- metadata$age - mean(metadata$age)

design_matrix <- model.matrix(
  ~ group + age_centered + sex + metformin + batch,
  data = metadata
)

if (qr(design_matrix)$rank != ncol(design_matrix)) {
  stop("The ANCOM-BC2 design matrix is not full rank.")
}

detection_prevalence <- rowMeans(count_matrix >= detection_threshold)
keep_taxon <- detection_prevalence >= minimum_detection_prevalence

filter_report <- data.frame(
  Genus = rownames(count_matrix),
  detection_threshold = detection_threshold,
  detection_prevalence = detection_prevalence,
  minimum_required_prevalence = minimum_detection_prevalence,
  retained = keep_taxon,
  stringsAsFactors = FALSE
)

if (sum(keep_taxon) < 2) {
  stop("The training filter retained fewer than two taxa.")
}

filtered_counts <- count_matrix[keep_taxon, , drop = FALSE]

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.table(
  filter_report,
  filter_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("======================================\n")
cat("ANCOM-BC2 TRAINING ANALYSIS\n")
cat("======================================\n\n")
cat("Input taxa:", nrow(count_matrix), "\n")
cat("Retained taxa:", nrow(filtered_counts), "\n")
cat("Excluded taxa:", paste(rownames(count_matrix)[!keep_taxon], collapse = ", "), "\n")
cat("Samples:", ncol(filtered_counts), "\n")
cat("Design-matrix rank:", qr(design_matrix)$rank, "of", ncol(design_matrix), "\n")
cat("Formula: group + age_centered + sex + metformin + batch\n\n")

set.seed(20260831)

fit <- ancombc2(
  data = filtered_counts,
  taxa_are_rows = TRUE,
  aggregate_data = filtered_counts,
  meta_data = metadata,
  fix_formula = "group + age_centered + sex + metformin + batch",
  rand_formula = NULL,
  p_adj_method = "BH",
  pseudo = 0,
  pseudo_sens = TRUE,
  prv_cut = 0,
  lib_cut = 0,
  s0_perc = 0.05,
  group = "group",
  struc_zero = TRUE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE
)

full_results <- fit$res

taxon_column <- intersect(c("taxon", "Genus"), names(full_results))
if (length(taxon_column) != 1) {
  stop(
    "Could not identify the taxon column in ANCOM-BC2 results. Columns: ",
    paste(names(full_results), collapse = ", ")
  )
}

names(full_results)[names(full_results) == taxon_column] <- "Genus"

target_suffix <- "groupT2D"
target_columns <- c(
  paste0("lfc_", target_suffix),
  paste0("se_", target_suffix),
  paste0("W_", target_suffix),
  paste0("p_", target_suffix),
  paste0("q_", target_suffix),
  paste0("diff_", target_suffix),
  paste0("passed_ss_", target_suffix),
  paste0("diff_robust_", target_suffix)
)

missing_result_columns <- setdiff(target_columns, names(full_results))
if (length(missing_result_columns) > 0) {
  stop(
    "Missing expected ANCOM-BC2 result columns: ",
    paste(missing_result_columns, collapse = ", "),
    ". Available columns: ", paste(names(full_results), collapse = ", ")
  )
}

group_results <- full_results[, c("Genus", target_columns), drop = FALSE]
names(group_results) <- c(
  "Genus",
  "estimated_log_effect",
  "standard_error",
  "W_statistic",
  "p_value",
  "q_value_BH",
  "differential",
  "passed_pseudocount_sensitivity",
  "differential_robust"
)

truth_group <- truth[, c("Genus", "group_T2D_log_effect")]
names(truth_group)[2] <- "simulated_true_log_effect"
group_results <- merge(group_results, truth_group, by = "Genus", all.x = TRUE)
group_results <- group_results[order(group_results$q_value_BH), ]

write.table(
  full_results,
  full_results_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  group_results,
  group_results_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nGroup-effect results ordered by BH q-value:\n")
print(group_results, row.names = FALSE, digits = 4)

cat("\nArtificial T2D signals encoded in the simulation:\n")
print(
  truth_group[truth_group$simulated_true_log_effect != 0, ],
  row.names = FALSE
)

cat("\nRobust discoveries for groupT2D:\n")
print(
  group_results$Genus[group_results$differential_robust %in% TRUE]
)

cat("\nFiles saved:\n")
cat("-", filter_path, "\n")
cat("-", full_results_path, "\n")
cat("-", group_results_path, "\n")

cat("\nTraining warning:\n")
cat("These simulated results are not evidence about T2D or method performance.\n")
