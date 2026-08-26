import numpy as np
import pandas as pd

# Load feature table and metadata
feature_table = pd.read_csv(
    "data/example_feature_table.tsv",
    sep="\t"
)

metadata = pd.read_csv(
    "data/example_metadata.tsv",
    sep="\t"
)

# Put Feature_ID in the index
feature_table = feature_table.set_index("Feature_ID")

results = []

# Calculate alpha-diversity indices for each sample
for sample_id in feature_table.columns:

    counts = feature_table[sample_id]

    # Remove absent features
    positive_counts = counts[counts > 0]

    # Observed richness
    observed_asvs = len(positive_counts)

    # Relative abundances among observed ASVs
    proportions = positive_counts / positive_counts.sum()

    # Shannon diversity
    shannon = -(proportions * np.log(proportions)).sum()

    # Simpson diversity expressed as 1-D
    simpson = 1 - (proportions ** 2).sum()

    results.append(
        {
            "sample_id": sample_id,
            "observed_asvs": observed_asvs,
            "shannon": shannon,
            "simpson_1_minus_D": simpson
        }
    )

# Convert results into a dataframe
alpha_diversity = pd.DataFrame(results)

# Add participant metadata
alpha_diversity = alpha_diversity.merge(
    metadata,
    on="sample_id",
    how="left",
    validate="one_to_one"
)

# Round indices for display
alpha_diversity["shannon"] = alpha_diversity["shannon"].round(3)
alpha_diversity["simpson_1_minus_D"] = (
    alpha_diversity["simpson_1_minus_D"].round(3)
)

print("======================================")
print("ALPHA DIVERSITY RESULTS")
print("======================================")
print()
print(alpha_diversity.to_string(index=False))

# Save results
output_path = "results/alpha_diversity.tsv"

alpha_diversity.to_csv(
    output_path,
    sep="\t",
    index=False
)

print()
print("Results saved to:")
print(output_path)
