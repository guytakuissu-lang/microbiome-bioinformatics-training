from pathlib import Path

import pandas as pd


INPUT_PATH = Path("results/genus_counts.tsv")
REPORT_PATH = Path("results/genus_prevalence_report.tsv")
FILTERED_PATH = Path("results/genus_counts_filtered_training.tsv")

# Pedagogical thresholds only; they are not universal recommendations.
DETECTION_THRESHOLD = 10
MINIMUM_PREVALENCE = 2 / 3


counts = pd.read_csv(INPUT_PATH, sep="\t", index_col=0)
counts = counts.apply(pd.to_numeric, errors="raise")

if counts.isna().any().any():
    raise ValueError("The genus count table contains missing values.")

if (counts < 0).any().any():
    raise ValueError("The genus count table contains negative values.")

number_of_samples = counts.shape[1]

if number_of_samples == 0:
    raise ValueError("The genus count table contains no samples.")

# A taxon is considered detected in a sample when its count is at least 10.
detected = counts >= DETECTION_THRESHOLD
samples_detected = detected.sum(axis=1)
prevalence = samples_detected / number_of_samples

report = pd.DataFrame(
    {
        "total_count": counts.sum(axis=1),
        "mean_count": counts.mean(axis=1),
        "samples_detected": samples_detected,
        "number_of_samples": number_of_samples,
        "prevalence_fraction": prevalence,
        "prevalence_percent": prevalence * 100,
        "keep_training_filter": prevalence >= MINIMUM_PREVALENCE,
    }
)

report.index.name = "Genus"

filtered_counts = counts.loc[report["keep_training_filter"]].copy()

if filtered_counts.empty:
    raise ValueError("The training filter removed every genus.")

REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
report.to_csv(REPORT_PATH, sep="\t", float_format="%.6f")
filtered_counts.to_csv(FILTERED_PATH, sep="\t")

retained_taxa = report.index[report["keep_training_filter"]].tolist()
excluded_taxa = report.index[~report["keep_training_filter"]].tolist()

print("======================================")
print("TAXON PREVALENCE AND TRAINING FILTER")
print("======================================")
print()
print(f"Detection threshold: count >= {DETECTION_THRESHOLD}")
print(
    "Minimum prevalence: "
    f"{MINIMUM_PREVALENCE:.3f} "
    f"({MINIMUM_PREVALENCE * 100:.1f}%)"
)
print(f"Number of samples: {number_of_samples}")

print()
print("Prevalence report:")
display_columns = [
    "total_count",
    "samples_detected",
    "prevalence_percent",
    "keep_training_filter",
]
print(report[display_columns].round(2).to_string())

print()
print(f"Retained genera ({len(retained_taxa)}):")
print(retained_taxa)

print()
print(f"Excluded genera ({len(excluded_taxa)}):")
print(excluded_taxa)

print()
print("Methodological warning:")
print("- These thresholds are arbitrary and used only for training.")
print("- A zero can reflect absence, undersampling, or technical failure.")
print("- Real filtering must be prespecified and method-appropriate.")
print("- Filtering can remove rare but biologically relevant taxa.")

print()
print(f"Report saved to: {REPORT_PATH}")
print(f"Filtered training table saved to: {FILTERED_PATH}")
