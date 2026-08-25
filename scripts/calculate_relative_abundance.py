import pandas as pd

# Load feature table
feature_table = pd.read_csv(
    "data/example_feature_table.tsv",
    sep="\t"
)

# Separate feature IDs from abundance data
feature_ids = feature_table["Feature_ID"]
abundance = feature_table.iloc[:, 1:]

# Calculate relative abundance
relative_abundance = abundance.div(
    abundance.sum(axis=0),
    axis=1
) * 100

# Restore feature IDs
relative_abundance.insert(
    0,
    "Feature_ID",
    feature_ids
)

print("======================================")
print("RELATIVE ABUNDANCE")
print("======================================")

print("\nRelative abundance (%):")
print(relative_abundance.round(2).to_string(index=False))

print("\nColumn sums (%):")
print(
    relative_abundance.iloc[:, 1:]
    .sum()
    .round(2)
)

print("\n======================================")
print("Analysis completed.")
print("======================================")
