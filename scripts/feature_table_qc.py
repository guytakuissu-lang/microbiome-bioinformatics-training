import pandas as pd

# Load feature table
feature_table = pd.read_csv(
    "data/example_feature_table.tsv",
    sep="\t"
)

print("======================================")
print("FEATURE TABLE QUALITY CONTROL")
print("======================================")

# Basic dimensions
print("\nFeature table dimensions:")
print(feature_table.shape)

# Feature IDs
print("\nFeature IDs:")
print(feature_table["Feature_ID"].tolist())

print("\nNumber of features:")
print(feature_table["Feature_ID"].nunique())

# Sample IDs
sample_ids = feature_table.columns[1:].tolist()

print("\nSample IDs:")
print(sample_ids)

print("\nNumber of samples:")
print(len(sample_ids))

# Missing values
print("\nMissing values:")
print(feature_table.isna().sum())

# Negative values
abundance = feature_table.iloc[:, 1:]

print("\nNegative values:")
print((abundance < 0).sum().sum())

# Duplicate feature IDs
print("\nDuplicate feature IDs:")
print(feature_table["Feature_ID"].duplicated().sum())

# Total counts per sample
print("\nTotal counts per sample:")
print(abundance.sum())

# Total counts per feature
print("\nTotal counts per feature:")
feature_totals = abundance.sum(axis=1)
feature_totals.index = feature_table["Feature_ID"]
print(feature_totals)

# Load metadata
metadata = pd.read_csv(
    "data/example_metadata.tsv",
    sep="\t"
)

metadata_ids = metadata["sample_id"].tolist()

print("\nMetadata sample IDs:")
print(metadata_ids)

print("\nSample IDs in feature table but not metadata:")
print(sorted(set(sample_ids) - set(metadata_ids)))

print("\nSample IDs in metadata but not feature table:")
print(sorted(set(metadata_ids) - set(sample_ids)))

print("\n======================================")
print("QC completed.")
print("======================================")
