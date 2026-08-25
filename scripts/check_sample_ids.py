import pandas as pd

metadata = pd.read_csv("data/example_metadata.tsv", sep="\t")

print("Sample ID quality control")
print("=========================")

print("Number of samples:", len(metadata))

print("All sample IDs unique:",
      metadata["sample_id"].is_unique)

print("Missing sample IDs:",
      metadata["sample_id"].isna().sum())

print("Empty sample IDs:",
      (metadata["sample_id"].str.strip() == "").sum())

print("\nSample IDs:")
print(metadata["sample_id"].tolist())
