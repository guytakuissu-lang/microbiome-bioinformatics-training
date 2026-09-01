suppressPackageStartupMessages(library(ggplot2))

recovery_path <- "results/ancombc2_training_signal_recovery.tsv"
excluded_path <- "results/ancombc2_training_excluded_true_signals.tsv"
output_path <- "figures/ancombc2_training_signal_recovery.png"

recovery <- read.delim(recovery_path, check.names = FALSE)
excluded <- read.delim(excluded_path, check.names = FALSE)

required_recovery_columns <- c(
  "Genus", "covariate", "estimated_log_effect",
  "approximate_ci_lower", "approximate_ci_upper",
  "simulated_true_log_effect", "simulated_true_signal",
  "differential_robust"
)
required_excluded_columns <- c(
  "Genus", "covariate", "simulated_true_log_effect", "retained"
)

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

check_columns(recovery, required_recovery_columns, "Recovery table")
check_columns(excluded, required_excluded_columns, "Excluded-signal table")

analyzed_signals <- recovery[recovery$simulated_true_signal %in% TRUE, ]
analyzed_signals$status <- ifelse(
  analyzed_signals$differential_robust %in% TRUE,
  "Robustly detected",
  "Analyzed but not robustly detected"
)

excluded_signals <- data.frame(
  Genus = excluded$Genus,
  covariate = excluded$covariate,
  estimated_log_effect = NA_real_,
  approximate_ci_lower = NA_real_,
  approximate_ci_upper = NA_real_,
  simulated_true_log_effect = excluded$simulated_true_log_effect,
  status = "Removed by prevalence filter",
  stringsAsFactors = FALSE
)

signal_data <- rbind(
  analyzed_signals[, names(excluded_signals), drop = FALSE],
  excluded_signals
)

if (nrow(signal_data) == 0) {
  stop("No simulated signals were available for plotting.")
}

covariate_order <- c(
  "T2D group", "Age", "Male sex", "Metformin", "Batch 2"
)
signal_data$covariate_order <- match(
  signal_data$covariate,
  covariate_order
)
signal_data <- signal_data[
  order(signal_data$covariate_order, signal_data$Genus),
]

signal_data$row_label <- paste(signal_data$covariate, signal_data$Genus, sep = " — ")
signal_data$row_label <- factor(
  signal_data$row_label,
  levels = rev(signal_data$row_label)
)

status_levels <- c(
  "Robustly detected",
  "Analyzed but not robustly detected",
  "Removed by prevalence filter"
)
signal_data$status <- factor(signal_data$status, levels = status_levels)

estimated_points <- signal_data[is.finite(signal_data$estimated_log_effect), ]

point_data <- rbind(
  data.frame(
    row_label = signal_data$row_label,
    log_effect = signal_data$simulated_true_log_effect,
    marker = "Simulated true effect",
    status = signal_data$status
  ),
  data.frame(
    row_label = estimated_points$row_label,
    log_effect = estimated_points$estimated_log_effect,
    marker = "ANCOM-BC2 estimate",
    status = estimated_points$status
  )
)
point_data$marker <- factor(
  point_data$marker,
  levels = c("ANCOM-BC2 estimate", "Simulated true effect")
)

status_colours <- c(
  "Robustly detected" = "#0072B2",
  "Analyzed but not robustly detected" = "#E69F00",
  "Removed by prevalence filter" = "#7A1F5C"
)

plot_object <- ggplot(signal_data, aes(y = row_label)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey45"
  ) +
  geom_segment(
    data = estimated_points,
    aes(
      x = approximate_ci_lower,
      xend = approximate_ci_upper,
      yend = row_label,
      colour = status
    ),
    linewidth = 1,
    alpha = 0.85
  ) +
  geom_point(
    data = point_data,
    aes(
      x = log_effect,
      y = row_label,
      colour = status,
      shape = marker
    ),
    size = 4.2,
    stroke = 1
  ) +
  scale_colour_manual(values = status_colours, name = "Outcome") +
  scale_shape_manual(
    values = c("ANCOM-BC2 estimate" = 16, "Simulated true effect" = 17),
    name = "Marker"
  ) +
  labs(
    title = "Three of six simulated signals were robustly detected",
    subtitle = paste0(
      "Two signals were analyzed but not robust; ",
      "one rare-taxon signal was removed before analysis"
    ),
    x = "Log effect in the covariate-specific unit",
    y = NULL,
    caption = paste(
      "Approximate intervals: estimate +/- 1.96 standard errors.",
      "Known simulation truth; 11 retained taxa;",
      "no method-performance or biological inference.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 12.5, margin = margin(b = 12)),
    plot.caption = element_text(
      size = 10.5,
      colour = "grey35",
      hjust = 0,
      margin = margin(t = 12)
    ),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(15, 20, 15, 15)
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_path,
  plot = plot_object,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Figure saved to:", output_path, "\n")
cat("\nSummary represented in the figure:\n")
print(table(signal_data$status, useNA = "ifany"))
cat("\nInterpretation warning:\n")
cat("- Detection status is defined using known simulation truth.\n")
cat("- Failure to detect a signal does not prove absence of an effect.\n")
cat("- Taxa removed by filtering cannot be evaluated by the fitted model.\n")
