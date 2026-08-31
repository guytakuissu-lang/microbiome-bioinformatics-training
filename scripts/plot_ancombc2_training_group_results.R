suppressPackageStartupMessages(library(ggplot2))

input_path <- "results/ancombc2_training_group_results.tsv"
table_path <- "results/ancombc2_training_group_interpretation.tsv"
figure_path <- "figures/ancombc2_training_group_forest.png"

results <- read.delim(input_path, check.names = FALSE)

required_columns <- c(
  "Genus", "estimated_log_effect", "standard_error", "q_value_BH",
  "passed_pseudocount_sensitivity", "differential_robust",
  "simulated_true_log_effect"
)

missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

results$ci_lower_log <- results$estimated_log_effect - 1.96 * results$standard_error
results$ci_upper_log <- results$estimated_log_effect + 1.96 * results$standard_error
results$estimated_abundance_ratio <- exp(results$estimated_log_effect)
results$ci_lower_abundance_ratio <- exp(results$ci_lower_log)
results$ci_upper_abundance_ratio <- exp(results$ci_upper_log)
results$simulated_true_abundance_ratio <- exp(results$simulated_true_log_effect)

results$result_status <- ifelse(
  results$differential_robust,
  "Robust BH q < 0.05",
  "Not robustly differential"
)

results <- results[order(results$estimated_log_effect), ]
results$Genus <- factor(results$Genus, levels = results$Genus)

dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(figure_path), recursive = TRUE, showWarnings = FALSE)

write.table(
  results,
  table_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_object <- ggplot(
  results,
  aes(y = Genus, x = estimated_abundance_ratio)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_errorbar(
    aes(
      xmin = ci_lower_abundance_ratio,
      xmax = ci_upper_abundance_ratio,
      color = result_status
    ),
    orientation = "y",
    width = 0.20,
    linewidth = 0.7
  ) +
  geom_point(
    aes(color = result_status),
    size = 2.8
  ) +
  geom_point(
    aes(x = simulated_true_abundance_ratio),
    shape = 17,
    size = 2.6,
    color = "black"
  ) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "Robust BH q < 0.05" = "#D55E00",
      "Not robustly differential" = "#0072B2"
    )
  ) +
  labs(
    title = "ANCOM-BC2 group effects in simulated training data",
    subtitle = paste(
      "Circles: estimates with approximate 95% intervals;",
      "triangles: simulated true effects"
    ),
    x = "Estimated abundance ratio, exp(log effect), on log scale",
    y = NULL,
    color = "Result",
    caption = paste(
      "T2D versus Control among metformin-unexposed participants,",
      "at mean age, Female, Batch1. Simulated data; no biological inference."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0, color = "grey35")
  )

ggsave(
  filename = figure_path,
  plot = plot_object,
  width = 9,
  height = 6.5,
  dpi = 300
)

cat("======================================\n")
cat("ANCOM-BC2 GROUP-EFFECT INTERPRETATION\n")
cat("======================================\n\n")

display_columns <- c(
  "Genus", "estimated_log_effect", "ci_lower_log", "ci_upper_log",
  "estimated_abundance_ratio", "ci_lower_abundance_ratio",
  "ci_upper_abundance_ratio", "q_value_BH", "differential_robust",
  "simulated_true_log_effect"
)

print(results[, display_columns], row.names = FALSE, digits = 4)

cat("\nCoefficient interpretation:\n")
cat(paste(
  "groupT2D compares T2D with Control among participants with metformin = No,",
  "at the mean age, sex = Female, and batch = Batch1.\n"
))

cat("\nFiles saved:\n")
cat("-", table_path, "\n")
cat("-", figure_path, "\n")

cat("\nTraining warning:\n")
cat("Approximate intervals and simulated signals do not validate ANCOM-BC2 performance.\n")
