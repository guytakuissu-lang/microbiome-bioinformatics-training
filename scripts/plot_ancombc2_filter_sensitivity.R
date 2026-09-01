suppressPackageStartupMessages(library(ggplot2))

input_path <- "results/ancombc2_training_filter_sensitivity_coefficients.tsv"
output_path <- "figures/ancombc2_training_filter_sensitivity.png"

results <- read.delim(input_path, check.names = FALSE)

required_columns <- c(
  "minimum_prevalence_percent", "Genus", "covariate", "retained",
  "q_value_BH", "differential_robust", "simulated_true_signal",
  "recovery_classification"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0) {
  stop(
    "Filter-sensitivity results are missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

signal_results <- results[results$simulated_true_signal %in% TRUE, ]
if (nrow(signal_results) == 0) {
  stop("No simulated signals were available for plotting.")
}

signal_results$plot_status <- ifelse(
  !signal_results$retained,
  "Removed by filter",
  ifelse(
    signal_results$differential_robust %in% TRUE,
    "Robustly detected",
    "Analyzed but not robustly detected"
  )
)

signal_results$cell_label <- ifelse(
  !signal_results$retained,
  "Filtered",
  paste0("q = ", formatC(signal_results$q_value_BH, format = "f", digits = 3))
)

signal_results$threshold_label <- paste0(
  formatC(signal_results$minimum_prevalence_percent, format = "f", digits = 0),
  "%"
)
signal_results$threshold_label <- factor(
  signal_results$threshold_label,
  levels = c("5%", "20%", "50%")
)

signal_results$row_label <- paste(
  signal_results$covariate,
  signal_results$Genus,
  sep = " — "
)
row_order <- c(
  "T2D group — Genus_03",
  "T2D group — Genus_07",
  "Age — Genus_04",
  "Male sex — Genus_05",
  "Metformin — Genus_10",
  "Batch 2 — Genus_12"
)

missing_rows <- setdiff(unique(signal_results$row_label), row_order)
if (length(missing_rows) > 0) {
  stop("Unexpected simulated-signal rows: ", paste(missing_rows, collapse = ", "))
}
signal_results$row_label <- factor(
  signal_results$row_label,
  levels = rev(row_order)
)

status_levels <- c(
  "Robustly detected",
  "Analyzed but not robustly detected",
  "Removed by filter"
)
signal_results$plot_status <- factor(
  signal_results$plot_status,
  levels = status_levels
)

status_colours <- c(
  "Robustly detected" = "#0072B2",
  "Analyzed but not robustly detected" = "#E69F00",
  "Removed by filter" = "#7A1F5C"
)

plot_object <- ggplot(
  signal_results,
  aes(x = threshold_label, y = row_label, fill = plot_status)
) +
  geom_tile(colour = "white", linewidth = 1.4, width = 0.96, height = 0.92) +
  geom_text(
    aes(label = cell_label),
    colour = "white",
    fontface = "bold",
    size = 4.3
  ) +
  scale_fill_manual(values = status_colours, name = "Signal status") +
  scale_x_discrete(position = "top") +
  labs(
    title = "Filtering changes analyzability, but not the three robust conclusions",
    subtitle = paste0(
      "Genus_12 is retained at 5% prevalence but its simulated batch effect remains non-robust"
    ),
    x = "Minimum prevalence threshold with detection defined as count >= 10",
    y = NULL,
    caption = paste(
      "Cell labels show BH-adjusted q-values; filtered taxa have no fitted estimate.",
      "Thresholds were prespecified for this training sensitivity analysis.",
      "Single simulation with 10-12 analyzed taxa: no biological or method-performance inference.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12.5, margin = margin(b = 14)),
    plot.caption = element_text(
      size = 10.5,
      colour = "grey35",
      hjust = 0,
      margin = margin(t = 14)
    ),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(margin = margin(b = 10)),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(15, 20, 15, 15)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_path,
  plot = plot_object,
  width = 11.5,
  height = 7.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Figure saved to:", output_path, "\n")
cat("\nSignal status counts across all thresholds:\n")
print(table(signal_results$plot_status, useNA = "ifany"))
cat("\nInterpretation warning:\n")
cat("- Retention does not guarantee robust detection.\n")
cat("- Removal by filtering is not statistical evidence of absence.\n")
cat("- Threshold stability in this simulation does not establish general robustness.\n")
