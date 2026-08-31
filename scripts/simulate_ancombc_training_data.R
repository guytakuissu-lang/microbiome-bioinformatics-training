set.seed(20260831)

output_dir <- "data"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_per_group <- 30L
n_samples <- 2L * n_per_group
taxa <- sprintf("Genus_%02d", 1:12)

sample_id <- sprintf("T%03d", seq_len(n_samples))
group <- factor(
  rep(c("Control", "T2D"), each = n_per_group),
  levels = c("Control", "T2D")
)

# Balance sex and sequencing batch within each disease group.
sex <- unlist(lapply(
  split(seq_len(n_samples), group),
  function(index) sample(rep(c("Female", "Male"), each = length(index) / 2))
))
sex <- factor(sex, levels = c("Female", "Male"))

batch <- unlist(lapply(
  split(seq_len(n_samples), group),
  function(index) sample(rep(c("Batch1", "Batch2"), each = length(index) / 2))
))
batch <- factor(batch, levels = c("Batch1", "Batch2"))

# Controls are unexposed; half of the simulated T2D participants receive
# metformin. This is a teaching design, not a recommendation for causal coding.
metformin <- rep("No", n_samples)
t2d_index <- which(group == "T2D")
metformin[sample(t2d_index, n_per_group / 2)] <- "Yes"
metformin <- factor(metformin, levels = c("No", "Yes"))

age <- numeric(n_samples)
age[group == "Control"] <- rnorm(n_per_group, mean = 52, sd = 7)
age[group == "T2D"] <- rnorm(n_per_group, mean = 56, sd = 7)
age <- pmin(pmax(round(age), 35), 75)

library_size <- sample(8000:20000, n_samples, replace = TRUE)

metadata <- data.frame(
  sample_id = sample_id,
  group = group,
  age = age,
  sex = sex,
  metformin = metformin,
  batch = batch,
  library_size = library_size,
  stringsAsFactors = FALSE
)

# Log-scale baseline abundance proxies. The final taxa are intentionally rare
# to create a realistic mixture of common and sparse features.
baseline_log_abundance <- c(
  2.4, 2.1, 1.8, 1.5, 1.2, 0.8,
  0.5, 0.2, -0.5, -1.0, -3.5, -5.0
)

# Known effects used only to evaluate the training workflow.
group_log_effect <- rep(0, length(taxa))
group_log_effect[3] <- 0.70
group_log_effect[7] <- -0.60

age_effect_per_10_years <- rep(0, length(taxa))
age_effect_per_10_years[4] <- 0.25

male_log_effect <- rep(0, length(taxa))
male_log_effect[5] <- 0.35

metformin_log_effect <- rep(0, length(taxa))
metformin_log_effect[10] <- 0.80

batch2_log_effect <- rep(0, length(taxa))
batch2_log_effect[12] <- 0.70

count_matrix <- matrix(
  0L,
  nrow = length(taxa),
  ncol = n_samples,
  dimnames = list(taxa, sample_id)
)

for (i in seq_len(n_samples)) {
  linear_predictor <- baseline_log_abundance

  if (group[i] == "T2D") {
    linear_predictor <- linear_predictor + group_log_effect
  }

  linear_predictor <- linear_predictor +
    age_effect_per_10_years * ((age[i] - 50) / 10)

  if (sex[i] == "Male") {
    linear_predictor <- linear_predictor + male_log_effect
  }

  if (metformin[i] == "Yes") {
    linear_predictor <- linear_predictor + metformin_log_effect
  }

  if (batch[i] == "Batch2") {
    linear_predictor <- linear_predictor + batch2_log_effect
  }

  # Sample-level biological heterogeneity on the latent abundance scale.
  linear_predictor <- linear_predictor + rnorm(length(taxa), mean = 0, sd = 0.45)

  abundance_proxy <- exp(linear_predictor)
  composition <- abundance_proxy / sum(abundance_proxy)

  count_matrix[, i] <- as.integer(
    rmultinom(1, size = library_size[i], prob = composition)
  )
}

count_table <- data.frame(
  Genus = rownames(count_matrix),
  count_matrix,
  check.names = FALSE
)

truth <- data.frame(
  Genus = taxa,
  group_T2D_log_effect = group_log_effect,
  age_per_10_years_log_effect = age_effect_per_10_years,
  sex_Male_log_effect = male_log_effect,
  metformin_Yes_log_effect = metformin_log_effect,
  batch_Batch2_log_effect = batch2_log_effect,
  stringsAsFactors = FALSE
)

# Verify that the planned multivariable design matrix is estimable.
design_matrix <- model.matrix(
  ~ group + age + sex + metformin + batch,
  data = metadata
)

if (qr(design_matrix)$rank != ncol(design_matrix)) {
  stop("The simulated design matrix is not full rank.")
}

if (!all(unname(colSums(count_matrix)) == as.numeric(library_size))) {
  stop("Count totals do not match the simulated library sizes.")
}

metadata_path <- file.path(output_dir, "simulated_ancombc_metadata.tsv")
counts_path <- file.path(output_dir, "simulated_ancombc_genus_counts.tsv")
truth_path <- file.path(output_dir, "simulated_ancombc_truth.tsv")

write.table(
  metadata,
  metadata_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  count_table,
  counts_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  truth,
  truth_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("======================================\n")
cat("ANCOM-BC2 TRAINING DATA SIMULATION\n")
cat("======================================\n\n")
cat("Samples:", n_samples, "\n")
cat("Taxa:", length(taxa), "\n")
cat("Design-matrix rank:", qr(design_matrix)$rank, "of", ncol(design_matrix), "\n")
cat("Library-size range:", min(library_size), "to", max(library_size), "\n")
cat("Total zero cells:", sum(count_matrix == 0), "of", length(count_matrix), "\n\n")

cat("Group counts:\n")
print(table(metadata$group))

cat("\nMetformin by group:\n")
print(table(metadata$group, metadata$metformin))

cat("\nKnown non-zero simulated effects:\n")
print(truth[rowSums(truth[, -1] != 0) > 0, ], row.names = FALSE)

cat("\nFiles saved:\n")
cat("-", metadata_path, "\n")
cat("-", counts_path, "\n")
cat("-", truth_path, "\n")

cat("\nTraining warning:\n")
cat("These data are simulated and must not be interpreted as evidence about T2D.\n")
