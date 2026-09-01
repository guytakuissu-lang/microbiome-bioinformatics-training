suppressPackageStartupMessages(library(ANCOMBC))

counts_path <- "data/simulated_ancombc_genus_counts.tsv"
metadata_path <- "data/simulated_ancombc_metadata.tsv"
truth_path <- "data/simulated_ancombc_truth.tsv"

taxon_filter_path <- "results/ancombc2_training_filter_sensitivity_taxa.tsv"
coefficient_path <- "results/ancombc2_training_filter_sensitivity_coefficients.tsv"
summary_path <- "results/ancombc2_training_filter_sensitivity_summary.tsv"

detection_threshold <- 10
minimum_prevalence_values <- c(0.05, 0.20, 0.50)

count_table <- read.delim(counts_path, check.names = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE)
truth <- read.delim(truth_path, check.names = FALSE)

required_metadata_columns <- c(
  "sample_id", "group", "age", "sex", "metformin", "batch"
)
required_truth_columns <- c(
  "Genus",
  "group_T2D_log_effect",
  "age_per_10_years_log_effect",
  "sex_Male_log_effect",
  "metformin_Yes_log_effect",
  "batch_Batch2_log_effect"
)

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

check_columns(count_table, "Genus", "Count table")
check_columns(metadata, required_metadata_columns, "Metadata")
check_columns(truth, required_truth_columns, "Truth table")

count_matrix <- as.matrix(count_table[, -1, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- count_table$Genus

if (!identical(colnames(count_matrix), metadata$sample_id)) {
  stop("Count-table columns and metadata sample IDs are not identically ordered.")
}
if (!setequal(rownames(count_matrix), truth$Genus)) {
  stop("Count-table genera and truth-table genera do not match.")
}

metadata$group <- factor(metadata$group, levels = c("Control", "T2D"))
metadata$sex <- factor(metadata$sex, levels = c("Female", "Male"))
metadata$metformin <- factor(metadata$metformin, levels = c("No", "Yes"))
metadata$batch <- factor(metadata$batch, levels = c("Batch1", "Batch2"))
metadata$age_centered <- metadata$age - mean(metadata$age)
rownames(metadata) <- metadata$sample_id

adjusted_design <- model.matrix(
  ~ group + age_centered + sex + metformin + batch,
  data = metadata
)
if (qr(adjusted_design)$rank != ncol(adjusted_design)) {
  stop("The adjusted design matrix is not full rank.")
}

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
  stringsAsFactors = FALSE
)

detection_prevalence <- rowMeans(count_matrix >= detection_threshold)

extract_threshold_results <- function(
  result_table,
  retained_taxa,
  minimum_prevalence
) {
  taxon_column <- intersect(c("taxon", "Genus"), names(result_table))
  if (length(taxon_column) != 1) {
    stop("Could not identify the taxon column in ANCOM-BC2 results.")
  }
  names(result_table)[names(result_table) == taxon_column] <- "Genus"

  coefficient_parts <- lapply(
    seq_len(nrow(coefficient_specification)),
    function(index) {
      specification_row <- coefficient_specification[index, , drop = FALSE]
      suffix <- specification_row$result_suffix
      scale_factor <- specification_row$estimate_scale
      expected_columns <- c(
        paste0("lfc_", suffix),
        paste0("se_", suffix),
        paste0("p_", suffix),
        paste0("q_", suffix),
        paste0("diff_robust_", suffix)
      )

      missing <- setdiff(expected_columns, names(result_table))
      if (length(missing) > 0) {
        stop(
          "Missing ANCOM-BC2 columns for ",
          specification_row$covariate,
          ": ",
          paste(missing, collapse = ", ")
        )
      }

      analyzed_truth_index <- match(result_table$Genus, truth$Genus)
      if (any(is.na(analyzed_truth_index))) {
        stop("At least one analyzed genus is absent from the truth table.")
      }

      analyzed <- data.frame(
        minimum_prevalence = minimum_prevalence,
        minimum_prevalence_percent = 100 * minimum_prevalence,
        Genus = result_table$Genus,
        covariate = specification_row$covariate,
        retained = TRUE,
        estimated_log_effect =
          result_table[[paste0("lfc_", suffix)]] * scale_factor,
        standard_error =
          result_table[[paste0("se_", suffix)]] * scale_factor,
        p_value = result_table[[paste0("p_", suffix)]],
        q_value_BH = result_table[[paste0("q_", suffix)]],
        differential_robust =
          result_table[[paste0("diff_robust_", suffix)]] %in% TRUE,
        simulated_true_log_effect = truth[[specification_row$truth_column]][
          analyzed_truth_index
        ],
        stringsAsFactors = FALSE
      )

      excluded_taxa <- setdiff(truth$Genus, retained_taxa)
      if (length(excluded_taxa) == 0) {
        excluded <- analyzed[0, , drop = FALSE]
      } else {
        excluded_truth_index <- match(excluded_taxa, truth$Genus)
        excluded <- data.frame(
          minimum_prevalence = rep(minimum_prevalence, length(excluded_taxa)),
          minimum_prevalence_percent = rep(
            100 * minimum_prevalence,
            length(excluded_taxa)
          ),
          Genus = excluded_taxa,
          covariate = rep(
            specification_row$covariate,
            length(excluded_taxa)
          ),
          retained = FALSE,
          estimated_log_effect = NA_real_,
          standard_error = NA_real_,
          p_value = NA_real_,
          q_value_BH = NA_real_,
          differential_robust = NA,
          simulated_true_log_effect = truth[[specification_row$truth_column]][
            excluded_truth_index
          ],
          stringsAsFactors = FALSE
        )
      }

      rbind(analyzed, excluded)
    }
  )

  threshold_results <- do.call(rbind, coefficient_parts)
  threshold_results$simulated_true_signal <-
    threshold_results$simulated_true_log_effect != 0
  threshold_results$recovery_classification <- ifelse(
    !threshold_results$retained,
    "Removed by filter",
    ifelse(
      threshold_results$simulated_true_signal &
        threshold_results$differential_robust %in% TRUE,
      "Detected simulated signal",
      ifelse(
        threshold_results$simulated_true_signal,
        "Missed simulated signal",
        ifelse(
          threshold_results$differential_robust %in% TRUE,
          "Robust discovery at simulated null",
          "Non-robust result at simulated null"
        )
      )
    )
  )
  threshold_results
}

taxon_filter_parts <- list()
coefficient_parts <- list()

cat("======================================\n")
cat("ANCOM-BC2 FILTER SENSITIVITY\n")
cat("======================================\n\n")
cat("Detection threshold: count >=", detection_threshold, "\n")
cat(
  "Minimum prevalence values:",
  paste0(100 * minimum_prevalence_values, "%", collapse = ", "),
  "\n"
)
cat("Adjusted formula: group + age_centered + sex + metformin + batch\n\n")

for (index in seq_along(minimum_prevalence_values)) {
  minimum_prevalence <- minimum_prevalence_values[index]
  retained <- detection_prevalence >= minimum_prevalence
  retained_taxa <- rownames(count_matrix)[retained]

  if (length(retained_taxa) < 2) {
    stop(
      "Fewer than two taxa were retained at minimum prevalence ",
      minimum_prevalence
    )
  }

  taxon_filter_parts[[index]] <- data.frame(
    minimum_prevalence = minimum_prevalence,
    minimum_prevalence_percent = 100 * minimum_prevalence,
    Genus = rownames(count_matrix),
    detection_threshold = detection_threshold,
    detection_prevalence = detection_prevalence,
    retained = retained,
    stringsAsFactors = FALSE
  )

  filtered_counts <- count_matrix[retained, , drop = FALSE]

  cat(
    "Running minimum prevalence",
    paste0(100 * minimum_prevalence, "%"),
    "with",
    nrow(filtered_counts),
    "taxa ...\n"
  )

  set.seed(20260901 + index)
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
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )

  coefficient_parts[[index]] <- extract_threshold_results(
    result_table = fit$res,
    retained_taxa = retained_taxa,
    minimum_prevalence = minimum_prevalence
  )
}

taxon_filter_results <- do.call(rbind, taxon_filter_parts)
coefficient_results <- do.call(rbind, coefficient_parts)

summarize_threshold <- function(minimum_prevalence) {
  subset_data <- coefficient_results[
    coefficient_results$minimum_prevalence == minimum_prevalence,
  ]
  taxon_subset <- taxon_filter_results[
    taxon_filter_results$minimum_prevalence == minimum_prevalence,
  ]

  data.frame(
    minimum_prevalence = minimum_prevalence,
    minimum_prevalence_percent = 100 * minimum_prevalence,
    retained_taxa = sum(taxon_subset$retained),
    simulated_signals_total = sum(subset_data$simulated_true_signal),
    simulated_signals_analyzed = sum(
      subset_data$retained & subset_data$simulated_true_signal
    ),
    simulated_signals_removed = sum(
      !subset_data$retained & subset_data$simulated_true_signal
    ),
    detected_simulated_signals = sum(
      subset_data$retained &
        subset_data$simulated_true_signal &
        subset_data$differential_robust %in% TRUE
    ),
    missed_simulated_signals = sum(
      subset_data$retained &
        subset_data$simulated_true_signal &
        !(subset_data$differential_robust %in% TRUE)
    ),
    robust_discoveries_at_simulated_null = sum(
      subset_data$retained &
        !subset_data$simulated_true_signal &
        subset_data$differential_robust %in% TRUE
    ),
    stringsAsFactors = FALSE
  )
}

sensitivity_summary <- do.call(
  rbind,
  lapply(minimum_prevalence_values, summarize_threshold)
)

coefficient_results <- coefficient_results[
  order(
    coefficient_results$minimum_prevalence,
    match(coefficient_results$covariate, coefficient_specification$covariate),
    coefficient_results$Genus
  ),
]

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.table(
  taxon_filter_results,
  taxon_filter_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  coefficient_results,
  coefficient_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  sensitivity_summary,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nSensitivity summary:\n")
print(sensitivity_summary, row.names = FALSE)

cat("\nGenus_12 batch-effect results:\n")
genus_12_batch <- coefficient_results[
  coefficient_results$Genus == "Genus_12" &
    coefficient_results$covariate == "Batch 2",
  c(
    "minimum_prevalence_percent", "retained",
    "estimated_log_effect", "q_value_BH", "differential_robust",
    "simulated_true_log_effect", "recovery_classification"
  )
]
print(genus_12_batch, row.names = FALSE, digits = 4)

cat("\nSimulated-signal status by threshold:\n")
signal_rows <- coefficient_results[
  coefficient_results$simulated_true_signal,
  c(
    "minimum_prevalence_percent", "Genus", "covariate", "retained",
    "estimated_log_effect", "q_value_BH", "differential_robust",
    "recovery_classification"
  )
]
print(signal_rows, row.names = FALSE, digits = 4)

cat("\nFiles saved:\n")
cat("-", taxon_filter_path, "\n")
cat("-", coefficient_path, "\n")
cat("-", summary_path, "\n")

cat("\nTraining warning:\n")
cat("- Thresholds were selected before inspecting these sensitivity results.\n")
cat("- A lower threshold can retain sparse, unstable taxa.\n")
cat("- A higher threshold can remove genuine simulated signals.\n")
cat("- This single small simulation does not estimate method performance.\n")
