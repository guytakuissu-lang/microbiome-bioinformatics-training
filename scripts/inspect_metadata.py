import pandas as pd

# Load metadata
metadata = pd.read_csv("data/example_metadata.tsv", sep="\t")

print("Metadata:")
print(metadata)

print("\nShape:")
print(metadata.shape)

print("\nColumns:")
print(metadata.columns.tolist())

print("\nData types:")
print(metadata.dtypes)

print("\nMissing values:")
print(metadata.isna().sum())

print("\nDuplicate sample IDs:")
print(metadata["sample_id"].duplicated().sum())

print("\nGroup counts:")
print(metadata["group"].value_counts())

print("\nSex counts:")
print(metadata["sex"].value_counts())

print("\nAge summary:")
print(metadata["age"].describe())
