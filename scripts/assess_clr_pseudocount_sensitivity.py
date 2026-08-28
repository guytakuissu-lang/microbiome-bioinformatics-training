from pathlib import Path

import numpy as np
import pandas as pd


INPUT_PATH = Path("results/genus_counts.tsv")
OUTPUT_PATH = Path("results/clr_pseudocount_sensitivity.tsv")
PSEUDOCOUNTS = [0.1, 0.5, 1.0]


counts = pd.read_csv(INPUT_PATH, sep="\t", index_col=0)
counts = counts.apply(pd.to_numeric, errors="raise")

if counts.isna().any().any():
    raise ValueError("The genus count table contains missing values.")

if (counts < 0).any().any():
    raise ValueError("The genus count table contains negative values.")

zero_locations = [
    (taxon, sample)
    for taxon in counts.index
    for sample in counts.columns
    if counts.loc[taxon, sample] == 0
]

all_results = []

for pseudocount in PSEUDOCOUNTS:
    adjusted_counts = counts + pseudocount
    log_counts = np.log(adjusted_counts)
    clr = log_counts.subtract(log_counts.mean(axis=0), axis=1)

    if not np.allclose(clr.sum(axis=0).to_numpy(), 0, atol=1e-10):
        raise ValueError(
            f"CLR coordinates do not sum to zero for {pseudocount}."
        )

    long_result = (
        clr.rename_axis("taxon")
        .reset_index()
        .melt(
            id_vars="taxon",
            var_name="sample_id",
            value_name="clr_value",
        )
    )
    long_result.insert(0, "pseudocount", pseudocount)
    all_results.append(long_result)

results = pd.concat(all_results, ignore_index=True)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(OUTPUT_PATH, sep="\t", index=False, float_format="%.6f")

print("======================================")
print("CLR PSEUDOCOUNT SENSITIVITY")
print("======================================")
print()
print(f"Pseudocounts evaluated: {PSEUDOCOUNTS}")
print(f"Observed zero locations: {zero_locations}")

if zero_locations:
    for taxon, sample_id in zero_locations:
        zero_summary = results.loc[
            (results["taxon"] == taxon)
            & (results["sample_id"] == sample_id),
            ["pseudocount", "clr_value"],
        ]

        print()
        print(f"Sensitivity for zero cell: {taxon}, {sample_id}")
        print(zero_summary.round(3).to_string(index=False))

        clr_range = (
            zero_summary["clr_value"].max()
            - zero_summary["clr_value"].min()
        )
        print(f"CLR range across pseudocounts: {clr_range:.3f}")

print()
print("Interpretation:")
print("- Different pseudocounts produce different CLR coordinates.")
print("- The largest sensitivity is expected for observed zeros.")
print("- No pseudocount is universally correct for real datasets.")
print("- Zero handling must match the downstream statistical method.")

print()
print(f"Results saved to: {OUTPUT_PATH}")
