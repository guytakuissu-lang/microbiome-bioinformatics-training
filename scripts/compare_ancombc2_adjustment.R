suppressPackageStartupMessages(library(ANCOMBC))

counts_path <- "data/simulated_ancombc_genus_counts.tsv"
metadata_path <- "data/simulated_ancombc_metadata.tsv"
truth_path <- "data/simulated_ancombc_truth.tsv"
filter_path <- "results/ancombc2_training_filter.tsv"
adjusted_results_path <- "results/ancombc2_training_group_results.tsv"
comparison_path <- "results/ancombc2_training_adjustment_comparison.tsv"

count_table <- read.delim(counts_path, check.names = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE)
truth <- read.delim(truth_path, check.names = FALSE)
filter_report <- read.delim(filter_path, check.names = FALSE)
adjusted_results <- read.delim(adjusted_results_path, check.names = FALSE)

required_count_columns <- "Genus"
required_metadata_columns <- c(
  "sample_id", "group", "age", "sex", "metformin", "batch"
)
required_truth_columns <- c(
  "Genus", "group_T2D_log_effect", "metformin_Yes_log_effect"
)
required_filter_columns <- c("Genus", "retained")
required_adjusted_columns <- c(
  "Genus", "estimated_log_effect", "standard_error", "p_value",
  "q_value_BH", "differential_robust"
)

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

check_columns(count_table, required_count_columns, "Count table")
check_columns(metadata, required_metadata_columns, "Metadata")
check_columns(truth, required_truth_columns, "Truth table")
check_columns(filter_report, required_filter_columns, "Filter report")
check_columns(
  adjusted_results,
  required_adjusted_columns,
  "Adjusted group-results table"
)

count_matrix <- as.matrix(count_table[, -1, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- count_table$Genus

if (!identical(colnames(count_matrix), metadata$sample_id)) {
  stop("Count-table columns and metadata sample IDs are not identically ordered.")
}

metadata$group <- factor(metadata$group, levels = c("Control", "T2D"))
metadata$sex <- factor(metadata$sex, levels = c("Female", "Male"))
metadata$metformin <- factor(metadata$metformin, levels = c("No", "Yes"))
metadata$batch <- factor(metadata$batch, levels = c("Batch1", "Batch2"))
metadata$age_centered <- metadata$age - mean(metadata$age)
rownames(metadata) <- metadata$sample_id

retained_taxa <- filter_report$Genus[filter_report$retained %in% TRUE]
missing_retained_taxa <- setdiff(retained_taxa, rownames(count_matrix))
if (length(missing_retained_taxa) > 0) {
  stop(
    "Retained taxa absent from the count table: ",
    paste(missing_retained_taxa, collapse = ", ")
  )
}

filtered_counts <- count_matrix[retained_taxa, , drop = FALSE]

unadjusted_design <- model.matrix(~ group, data = metadata)
adjusted_design <- model.matrix(
  ~ group + age_centered + sex + metformin + batch,
  data = metadata
)

if (qr(unadjusted_design)$rank != ncol(unadjusted_design)) {
  stop("The unadjusted design matrix is not full rank.")
}
if (qr(adjusted_design)$rank != ncol(adjusted_design)) {
  stop("The adjusted design matrix is not full rank.")
}

cat("======================================\n")
cat("ANCOM-BC2 ADJUSTMENT COMPARISON\n")
cat("======================================\n\n")
cat("Samples:", ncol(filtered_counts), "\n")
cat("Retained taxa:", nrow(filtered_counts), "\n")
cat("Unadjusted formula: group\n")
cat("Adjusted formula: group + age_centered + sex + metformin + batch\n")
cat(
  "Important overlap limitation: metformin exposure occurs only in the T2D group.\n\n"
)

set.seed(20260831)

unadjusted_fit <- ancombc2(
  data = filtered_counts,
  taxa_are_rows = TRUE,
  aggregate_data = filtered_counts,
  meta_data = metadata,
  fix_formula = "group",
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

extract_group_result <- function(result_table) {
  taxon_column <- intersect(c("taxon", "Genus"), names(result_table))
  if (length(taxon_column) != 1) {
    stop(
      "Could not identify the taxon column. Available columns: ",
      paste(names(result_table), collapse = ", ")
    )
  }

  names(result_table)[names(result_table) == taxon_column] <- "Genus"
  target_suffix <- "groupT2D"
  target_columns <- c(
    paste0("lfc_", target_suffix),
    paste0("se_", target_suffix),
    paste0("p_", target_suffix),
    paste0("q_", target_suffix),
    paste0("diff_robust_", target_suffix)
  )

  missing <- setdiff(target_columns, names(result_table))
  if (length(missing) > 0) {
    stop(
      "ANCOM-BC2 results are missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  extracted <- result_table[, c("Genus", target_columns), drop = FALSE]
  names(extracted) <- c(
    "Genus",
    "estimated_log_effect",
    "standard_error",
    "p_value",
    "q_value_BH",
    "differential_robust"
  )
  extracted
}

unadjusted_results <- extract_group_result(unadjusted_fit$res)

names(unadjusted_results)[-1] <- paste0(
  "unadjusted_", names(unadjusted_results)[-1]
)
adjusted_results <- adjusted_results[, required_adjusted_columns, drop = FALSE]
names(adjusted_results)[-1] <- paste0(
  "adjusted_", names(adjusted_results)[-1]
)

comparison <- merge(
  unadjusted_results,
  adjusted_results,
  by = "Genus",
  all = TRUE,
  sort = FALSE
)

truth_subset <- truth[, required_truth_columns, drop = FALSE]
names(truth_subset) <- c(
  "Genus",
  "simulated_group_T2D_log_effect",
  "simulated_metformin_Yes_log_effect"
)
comparison <- merge(comparison, truth_subset, by = "Genus", all.x = TRUE)

comparison$estimate_change_after_adjustment <-
  comparison$adjusted_estimated_log_effect -
  comparison$unadjusted_estimated_log_effect
comparison$absolute_estimate_change <- abs(
  comparison$estimate_change_after_adjustment
)
comparison$unadjusted_only_robust <-
  comparison$unadjusted_differential_robust %in% TRUE &
  !(comparison$adjusted_differential_robust %in% TRUE)
comparison$adjusted_only_robust <-
  comparison$adjusted_differential_robust %in% TRUE &
  !(comparison$unadjusted_differential_robust %in% TRUE)
comparison$simulated_metformin_only_signal <-
  comparison$simulated_group_T2D_log_effect == 0 &
  comparison$simulated_metformin_Yes_log_effect != 0

comparison <- comparison[
  order(comparison$absolute_estimate_change, decreasing = TRUE),
]

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.table(
  comparison,
  comparison_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

display_columns <- c(
  "Genus",
  "unadjusted_estimated_log_effect",
  "unadjusted_q_value_BH",
  "unadjusted_differential_robust",
  "adjusted_estimated_log_effect",
  "adjusted_q_value_BH",
  "adjusted_differential_robust",
  "estimate_change_after_adjustment",
  "simulated_group_T2D_log_effect",
  "simulated_metformin_Yes_log_effect"
)

cat("Comparison ordered by absolute estimate change:\n")
print(comparison[, display_columns], row.names = FALSE, digits = 4)

cat("\nGenus_10 pedagogical control:\n")
print(
  comparison[comparison$Genus == "Genus_10", display_columns],
  row.names = FALSE,
  digits = 4
)

cat("\nRobust discoveries in the unadjusted model:\n")
print(comparison$Genus[comparison$unadjusted_differential_robust %in% TRUE])

cat("\nRobust discoveries in the adjusted model:\n")
print(comparison$Genus[comparison$adjusted_differential_robust %in% TRUE])

cat("\nResults saved to:", comparison_path, "\n")
cat("\nInterpretation warning:\n")
cat("- Differences between models illustrate adjustment, not causality.\n")
cat("- Metformin has no exposed Control participants in this simulation.\n")
cat("- The adjusted group coefficient is therefore a conditional contrast.\n")
cat("- Simulated results do not validate ANCOM-BC2 performance.\n")
