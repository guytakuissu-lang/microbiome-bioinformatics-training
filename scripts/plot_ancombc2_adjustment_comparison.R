suppressPackageStartupMessages(library(ggplot2))

input_path <- "results/ancombc2_training_adjustment_comparison.tsv"
output_path <- "figures/ancombc2_training_adjustment_comparison.png"

comparison <- read.delim(input_path, check.names = FALSE)

required_columns <- c(
  "Genus",
  "unadjusted_estimated_log_effect",
  "adjusted_estimated_log_effect",
  "absolute_estimate_change",
  "simulated_group_T2D_log_effect",
  "simulated_metformin_Yes_log_effect"
)

missing_columns <- setdiff(required_columns, names(comparison))
if (length(missing_columns) > 0) {
  stop(
    "Comparison table is missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (any(!is.finite(comparison$unadjusted_estimated_log_effect)) ||
    any(!is.finite(comparison$adjusted_estimated_log_effect))) {
  stop("At least one model estimate is missing or non-finite.")
}

comparison <- comparison[
  order(comparison$absolute_estimate_change, decreasing = TRUE),
]

comparison$signal_type <- "No simulated group/metformin signal"
comparison$signal_type[
  comparison$simulated_group_T2D_log_effect != 0
] <- "Simulated group signal"
comparison$signal_type[
  comparison$simulated_group_T2D_log_effect == 0 &
    comparison$simulated_metformin_Yes_log_effect != 0
] <- "Simulated metformin-only signal"

comparison$signal_type <- factor(
  comparison$signal_type,
  levels = c(
    "No simulated group/metformin signal",
    "Simulated group signal",
    "Simulated metformin-only signal"
  )
)

comparison$Genus <- factor(
  comparison$Genus,
  levels = rev(as.character(comparison$Genus))
)

comparison$unadjusted_abundance_ratio <- exp(
  comparison$unadjusted_estimated_log_effect
)
comparison$adjusted_abundance_ratio <- exp(
  comparison$adjusted_estimated_log_effect
)

long_data <- rbind(
  data.frame(
    Genus = comparison$Genus,
    model = "Unadjusted",
    abundance_ratio = comparison$unadjusted_abundance_ratio,
    signal_type = comparison$signal_type
  ),
  data.frame(
    Genus = comparison$Genus,
    model = "Adjusted",
    abundance_ratio = comparison$adjusted_abundance_ratio,
    signal_type = comparison$signal_type
  )
)

long_data$model <- factor(
  long_data$model,
  levels = c("Unadjusted", "Adjusted")
)

signal_colours <- c(
  "No simulated group/metformin signal" = "#9E9E9E",
  "Simulated group signal" = "#0072B2",
  "Simulated metformin-only signal" = "#D55E00"
)

plot_object <- ggplot(comparison, aes(y = Genus)) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey45"
  ) +
  geom_segment(
    aes(
      x = unadjusted_abundance_ratio,
      xend = adjusted_abundance_ratio,
      yend = Genus,
      colour = signal_type
    ),
    linewidth = 1.1,
    alpha = 0.85
  ) +
  geom_point(
    data = long_data,
    aes(
      x = abundance_ratio,
      y = Genus,
      shape = model,
      fill = signal_type
    ),
    colour = "black",
    size = 4.2,
    stroke = 0.7
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = c("0.50", "0.75", "1.00", "1.50", "2.00")
  ) +
  scale_colour_manual(values = signal_colours, name = "Simulated signal") +
  scale_fill_manual(values = signal_colours, name = "Simulated signal") +
  scale_shape_manual(
    values = c("Unadjusted" = 21, "Adjusted" = 24),
    name = "Model"
  ) +
  labs(
    title = "Adjustment changes the estimated T2D association most for Genus_10",
    subtitle = paste0(
      "Unadjusted: group only; adjusted: group + age + sex + metformin + batch"
    ),
    x = "Estimated abundance ratio, exp(log effect), on log scale",
    y = NULL,
    caption = paste0(
      "Simulated data: no biological or causal inference. ",
      "Metformin exposure occurs only among T2D participants."
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
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = output_path,
  plot = plot_object,
  width = 12,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Figure saved to:", output_path, "\n")
cat("\nInterpretation:\n")
cat("- Horizontal segments connect unadjusted and adjusted estimates.\n")
cat("- Longer segments indicate greater sensitivity to covariate adjustment.\n")
cat("- Genus_10 is a simulated metformin-only control, not a biological result.\n")
cat("- Model differences do not establish causality or identify a confounder by themselves.\n")
