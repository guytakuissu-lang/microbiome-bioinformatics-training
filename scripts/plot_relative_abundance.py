import pandas as pd
import matplotlib.pyplot as plt

# Load feature table
feature_table = pd.read_csv(
    "data/example_feature_table.tsv",
    sep="\t"
)

# Load metadata
metadata = pd.read_csv(
    "data/example_metadata.tsv",
    sep="\t"
)

# Separate feature IDs from abundance data
feature_ids = feature_table["Feature_ID"]
abundance = feature_table.iloc[:, 1:].copy()

# Calculate relative abundance (%)
relative_abundance = abundance.div(
    abundance.sum(axis=0),
    axis=1
) * 100

# Use sample IDs as columns
relative_abundance.index = feature_ids

# Reorder samples according to metadata
sample_order = metadata["sample_id"].tolist()
relative_abundance = relative_abundance[sample_order]

# Create plot
fig, ax = plt.subplots(figsize=(10, 6))

bottom = pd.Series(0.0, index=sample_order)

for feature in relative_abundance.index:
    values = relative_abundance.loc[feature]

    ax.bar(
        sample_order,
        values,
        bottom=bottom,
        label=feature
    )

    bottom += values

# Add labels
ax.set_xlabel("Sample")
ax.set_ylabel("Relative abundance (%)")
ax.set_title("Microbial community composition")
ax.set_ylim(0, 100)

# Add group labels below samples
groups = metadata.set_index("sample_id").loc[sample_order, "group"]

for i, sample in enumerate(sample_order):
    ax.text(
        i,
        -7,
        groups[sample],
        ha="center",
        va="top"
    )

ax.legend(
    title="ASV",
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()

# Save figure
plt.savefig(
    "results/relative_abundance.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Relative abundance plot saved to:")
print("results/relative_abundance.png")
