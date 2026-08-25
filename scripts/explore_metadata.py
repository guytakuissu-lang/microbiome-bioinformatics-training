import pandas as pd

# Load metadata
metadata = pd.read_csv("data/example_metadata.tsv", sep="\t")

print("======================================")
print("METADATA EXPLORATORY ANALYSIS")
print("======================================")

# Basic information
print("\nDataset dimensions:")
print(metadata.shape)

print("\nVariables:")
print(metadata.columns.tolist())

# Group distribution
print("\nGroup distribution:")
print(metadata["group"].value_counts())

print("\nGroup proportions:")
print(metadata["group"].value_counts(normalize=True).round(3))

# Sex distribution
print("\nSex distribution:")
print(metadata["sex"].value_counts())

print("\nSex proportions:")
print(metadata["sex"].value_counts(normalize=True).round(3))

# Age statistics
print("\nOverall age statistics:")
print(metadata["age"].describe())

# Age by group
print("\nAge by group:")
print(
    metadata.groupby("group")["age"]
    .agg(["count", "mean", "std", "min", "max"])
)

# Sex by group
print("\nSex by group:")
print(
    pd.crosstab(metadata["group"], metadata["sex"])
)

print("\n======================================")
print("Analysis completed.")
print("======================================")
