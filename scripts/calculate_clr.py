from pathlib import Path

import numpy as np
import pandas as pd


INPUT_PATH = Path("results/genus_counts.tsv")
OUTPUT_PATH = Path("results/genus_clr_pseudocount_0.5.tsv")
PSEUDOCOUNT = 0.5


counts = pd.read_csv(INPUT_PATH, sep="\t", index_col=0)
counts = counts.apply(pd.to_numeric, errors="raise")

if counts.isna().any().any():
    raise ValueError("The genus count table contains missing values.")

if (counts < 0).any().any():
    raise ValueError("The genus count table contains negative values.")

if (counts.sum(axis=0) == 0).any():
    zero_samples = counts.columns[counts.sum(axis=0) == 0].tolist()
    raise ValueError(f"Samples with zero total counts: {zero_samples}")

# Locate observed zeros before replacement.
zero_locations = [
    (taxon, sample)
    for taxon in counts.index
    for sample in counts.columns
    if counts.loc[taxon, sample] == 0
]

# CLR cannot be calculated when values are zero because log(0) is undefined.
# Adding 0.5 to every cell is used here only as a transparent training example.
adjusted_counts = counts + PSEUDOCOUNT
log_counts = np.log(adjusted_counts)

# For each sample, subtract the mean of the log abundances. This is equivalent
# to taking the log of each abundance divided by the sample geometric mean.
clr = log_counts.subtract(log_counts.mean(axis=0), axis=1)

# A correctly calculated CLR composition sums to approximately zero per sample.
clr_sums = clr.sum(axis=0)

if not np.allclose(clr_sums.to_numpy(), 0, atol=1e-10):
    raise ValueError("CLR coordinates do not sum to zero per sample.")

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
clr.to_csv(OUTPUT_PATH, sep="\t", float_format="%.6f")

print("======================================")
print("CENTERED LOG-RATIO TRANSFORMATION")
print("======================================")
print()
print(f"Input: {INPUT_PATH}")
print(f"Pseudocount added to every cell: {PSEUDOCOUNT}")
print(f"Number of observed zeros: {len(zero_locations)}")

if zero_locations:
    print("Zero locations:")
    for taxon, sample in zero_locations:
        print(f"- {taxon}, {sample}")

print()
print("CLR coordinates:")
print(clr.round(3).to_string())

print()
print("CLR column sums (expected approximately zero):")
print(clr_sums.to_string())

print()
print("Interpretation reminder:")
print("- CLR values are log-ratio coordinates, not percentages.")
print("- A positive value is above the sample geometric mean.")
print("- A negative value is below the sample geometric mean.")
print("- Results involving zeros can depend strongly on zero handling.")

print()
print(f"Results saved to: {OUTPUT_PATH}")
