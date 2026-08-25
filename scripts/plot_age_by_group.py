import pandas as pd
import matplotlib.pyplot as plt

# Load metadata
metadata = pd.read_csv("data/example_metadata.tsv", sep="\t")

# Create figure
plt.figure(figsize=(7, 5))

# Plot age for each participant
for group in metadata["group"].unique():
    subset = metadata[metadata["group"] == group]
    plt.scatter(
        [group] * len(subset),
        subset["age"],
        s=80,
        label=group
    )

# Labels and title
plt.xlabel("Group")
plt.ylabel("Age (years)")
plt.title("Age distribution by group")

# Save figure
plt.tight_layout()
plt.savefig("results/age_by_group.png", dpi=300)

# Display figure
