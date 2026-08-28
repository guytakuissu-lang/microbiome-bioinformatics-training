from pathlib import Path

import pandas as pd


INPUT_PATH = Path("results/genus_counts.tsv")
OUTPUT_PATH = Path("results/genus_filter_sensitivity.tsv")

DETECTION_THRESHOLDS = [1, 5, 10, 20]
MINIMUM_PREVALENCES = [1 / 3, 2 / 3, 1.0]


counts = pd.read_csv(INPUT_PATH, sep="\t", index_col=0)
counts = counts.apply(pd.to_numeric, errors="raise")

if counts.isna().any().any():
    raise ValueError("The genus count table contains missing values.")

if (counts < 0).any().any():
    raise ValueError("The genus count table contains negative values.")

number_of_samples = counts.shape[1]

if number_of_samples == 0:
    raise ValueError("The genus count table contains no samples.")

records = []

for detection_threshold in DETECTION_THRESHOLDS:
    detected = counts >= detection_threshold
    prevalence = detected.mean(axis=1)

    for minimum_prevalence in MINIMUM_PREVALENCES:
        keep = prevalence >= minimum_prevalence
        retained = counts.index[keep].tolist()
        excluded = counts.index[~keep].tolist()

        records.append(
            {
                "detection_threshold": detection_threshold,
                "minimum_prevalence_fraction": minimum_prevalence,
                "minimum_prevalence_percent": minimum_prevalence * 100,
                "number_retained": len(retained),
                "number_excluded": len(excluded),
                "retained_genera": ";".join(retained),
                "excluded_genera": ";".join(excluded),
            }
        )

sensitivity = pd.DataFrame(records)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
sensitivity.to_csv(
    OUTPUT_PATH,
    sep="\t",
    index=False,
    float_format="%.6f",
)

print("======================================")
print("FILTER THRESHOLD SENSITIVITY")
print("======================================")
print()
print(f"Number of samples: {number_of_samples}")
print(f"Detection thresholds: {DETECTION_THRESHOLDS}")
print(
    "Minimum prevalence thresholds: "
    f"{[round(value * 100, 1) for value in MINIMUM_PREVALENCES]}%"
)

print()
print("Number of retained genera:")
summary = sensitivity.pivot(
    index="detection_threshold",
    columns="minimum_prevalence_percent",
    values="number_retained",
)
summary.columns = [f"{column:.1f}%" for column in summary.columns]
print(summary.to_string())

print()
print("Excluded genera for each threshold combination:")
for row in sensitivity.itertuples(index=False):
    excluded = row.excluded_genera if row.excluded_genera else "None"
    print(
        f"count >= {row.detection_threshold}; "
        f"prevalence >= {row.minimum_prevalence_percent:.1f}%: "
        f"{excluded}"
    )

print()
print("Interpretation:")
print("- The retained feature set changes when thresholds change.")
print("- Small sample sizes produce coarse prevalence increments.")
print("- Filtering choices must be justified and reported.")
print("- Sensitivity analyses should accompany consequential filters.")

print()
print(f"Results saved to: {OUTPUT_PATH}")
