full_results_path <- "results/ancombc2_training_full_results.tsv"
truth_path <- "data/simulated_ancombc_truth.tsv"
filter_path <- "results/ancombc2_training_filter.tsv"

recovery_path <- "results/ancombc2_training_signal_recovery.tsv"
summary_path <- "results/ancombc2_training_signal_recovery_summary.tsv"
excluded_signals_path <- "results/ancombc2_training_excluded_true_signals.tsv"

full_results <- read.delim(full_results_path, check.names = FALSE)
truth <- read.delim(truth_path, check.names = FALSE)
filter_report <- read.delim(filter_path, check.names = FALSE)

required_truth_columns <- c(
  "Genus",
  "group_T2D_log_effect",
  "age_per_10_years_log_effect",
  "sex_Male_log_effect",
  "metformin_Yes_log_effect",
  "batch_Batch2_log_effect"
)
required_filter_columns <- c("Genus", "retained")

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

check_columns(truth, required_truth_columns, "Truth table")
check_columns(filter_report, required_filter_columns, "Filter report")

taxon_column <- intersect(c("taxon", "Genus"), names(full_results))
if (length(taxon_column) != 1) {
  stop(
    "Could not identify the taxon column in the full ANCOM-BC2 results."
  )
}
names(full_results)[names(full_results) == taxon_column] <- "Genus"

coefficient_specification <- data.frame(
  covariate = c("T2D group", "Age", "Male sex", "Metformin", "Batch 2"),
  result_suffix = c(
    "groupT2D",
    "age_centered",
    "sexMale",
    "metforminYes",
    "batchBatch2"
  ),
  truth_column = c(
    "group_T2D_log_effect",
    "age_per_10_years_log_effect",
    "sex_Male_log_effect",
    "metformin_Yes_log_effect",
    "batch_Batch2_log_effect"
  ),
  estimate_scale = c(1, 10, 1, 1, 1),
  effect_unit = c(
    "T2D versus Control",
    "per 10-year age increase",
    "Male versus Female",
    "Yes versus No",
    "Batch2 versus Batch1"
  ),
  stringsAsFactors = FALSE
)

extract_coefficient <- function(specification_row) {
  suffix <- specification_row$result_suffix
  scale_factor <- specification_row$estimate_scale
  expected_columns <- c(
    paste0("lfc_", suffix),
    paste0("se_", suffix),
    paste0("p_", suffix),
    paste0("q_", suffix),
    paste0("diff_robust_", suffix)
  )

  missing <- setdiff(expected_columns, names(full_results))
  if (length(missing) > 0) {
    stop(
      "Full ANCOM-BC2 results are missing columns for ",
      specification_row$covariate,
      ": ",
      paste(missing, collapse = ", ")
    )
  }

  truth_index <- match(full_results$Genus, truth$Genus)
  if (any(is.na(truth_index))) {
    stop("At least one analyzed genus is absent from the truth table.")
  }

  estimated <- full_results[[paste0("lfc_", suffix)]] * scale_factor
  standard_error <- full_results[[paste0("se_", suffix)]] * scale_factor
  simulated_truth <- truth[[specification_row$truth_column]][truth_index]
  robust <- full_results[[paste0("diff_robust_", suffix)]] %in% TRUE
  true_signal <- simulated_truth != 0

  classification <- ifelse(
    true_signal & robust,
    "Detected simulated signal",
    ifelse(
      true_signal & !robust,
      "Missed simulated signal",
      ifelse(
        !true_signal & robust,
        "Robust discovery at simulated null",
        "Non-robust result at simulated null"
      )
    )
  )

  data.frame(
    Genus = full_results$Genus,
    covariate = specification_row$covariate,
    effect_unit = specification_row$effect_unit,
    estimated_log_effect = estimated,
    standard_error = standard_error,
    approximate_ci_lower = estimated - 1.96 * standard_error,
    approximate_ci_upper = estimated + 1.96 * standard_error,
    p_value = full_results[[paste0("p_", suffix)]],
    q_value_BH = full_results[[paste0("q_", suffix)]],
    differential_robust = robust,
    simulated_true_log_effect = simulated_truth,
    simulated_true_signal = true_signal,
    estimate_error = estimated - simulated_truth,
    absolute_estimate_error = abs(estimated - simulated_truth),
    sign_matches_truth = ifelse(
      true_signal,
      sign(estimated) == sign(simulated_truth),
      NA
    ),
    truth_within_approximate_interval =
      simulated_truth >= estimated - 1.96 * standard_error &
      simulated_truth <= estimated + 1.96 * standard_error,
    recovery_classification = classification,
    stringsAsFactors = FALSE
  )
}

recovery_parts <- lapply(
  seq_len(nrow(coefficient_specification)),
  function(index) {
    extract_coefficient(coefficient_specification[index, , drop = FALSE])
  }
)
recovery <- do.call(rbind, recovery_parts)

retained_taxa <- filter_report$Genus[filter_report$retained %in% TRUE]
if (!setequal(unique(recovery$Genus), retained_taxa)) {
  stop("Analyzed genera do not match the genera marked as retained.")
}

summarize_covariate <- function(covariate_name) {
  subset_data <- recovery[recovery$covariate == covariate_name, ]
  data.frame(
    covariate = covariate_name,
    analyzed_taxa = nrow(subset_data),
    simulated_signals = sum(subset_data$simulated_true_signal),
    robust_discoveries = sum(subset_data$differential_robust),
    detected_simulated_signals = sum(
      subset_data$simulated_true_signal & subset_data$differential_robust
    ),
    missed_simulated_signals = sum(
      subset_data$simulated_true_signal & !subset_data$differential_robust
    ),
    robust_discoveries_at_simulated_null = sum(
      !subset_data$simulated_true_signal & subset_data$differential_robust
    ),
    median_absolute_estimate_error = median(
      subset_data$absolute_estimate_error,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

recovery_summary <- do.call(
  rbind,
  lapply(coefficient_specification$covariate, summarize_covariate)
)

truth_long_parts <- lapply(
  seq_len(nrow(coefficient_specification)),
  function(index) {
    specification_row <- coefficient_specification[index, , drop = FALSE]
    data.frame(
      Genus = truth$Genus,
      covariate = specification_row$covariate,
      effect_unit = specification_row$effect_unit,
      simulated_true_log_effect = truth[[specification_row$truth_column]],
      stringsAsFactors = FALSE
    )
  }
)
truth_long <- do.call(rbind, truth_long_parts)
truth_long$retained <- filter_report$retained[
  match(truth_long$Genus, filter_report$Genus)
]

if (any(is.na(truth_long$retained))) {
  stop("At least one truth-table genus is absent from the filter report.")
}

excluded_true_signals <- truth_long[
  !(truth_long$retained %in% TRUE) &
    truth_long$simulated_true_log_effect != 0,
]

recovery <- recovery[
  order(
    match(recovery$covariate, coefficient_specification$covariate),
    !recovery$simulated_true_signal,
    recovery$q_value_BH
  ),
]

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.table(
  recovery,
  recovery_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  recovery_summary,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  excluded_true_signals,
  excluded_signals_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("======================================\n")
cat("ANCOM-BC2 SIMULATED-SIGNAL RECOVERY\n")
cat("======================================\n\n")

cat("Recovery summary by covariate:\n")
print(recovery_summary, row.names = FALSE, digits = 4)

cat("\nAnalyzed simulated signals and robust discoveries at simulated nulls:\n")
display_rows <- recovery[
  recovery$simulated_true_signal | recovery$differential_robust,
  c(
    "Genus", "covariate", "estimated_log_effect",
    "simulated_true_log_effect", "q_value_BH",
    "differential_robust", "recovery_classification"
  )
]
print(display_rows, row.names = FALSE, digits = 4)

cat("\nSimulated signals removed by the prevalence filter:\n")
if (nrow(excluded_true_signals) == 0) {
  cat("None\n")
} else {
  print(excluded_true_signals, row.names = FALSE, digits = 4)
}

cat("\nFiles saved:\n")
cat("-", recovery_path, "\n")
cat("-", summary_path, "\n")
cat("-", excluded_signals_path, "\n")

cat("\nTraining warning:\n")
cat("- Recovery labels use known simulation truth and are not available in real data.\n")
cat("- Approximate intervals are estimate +/- 1.96 standard errors.\n")
cat("- Eleven analyzed taxa are insufficient for method-performance validation.\n")
cat("- Filtering can remove a genuine simulated signal, as illustrated here.\n")
