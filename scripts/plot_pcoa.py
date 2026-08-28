from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


INPUT_FILE = Path("results/pcoa_coordinates.tsv")
OUTPUT_FILE = Path("figures/pcoa_bray_curtis.png")


# Load the coordinates produced by calculate_beta_diversity.py
pcoa = pd.read_csv(INPUT_FILE, sep="\t")

required_columns = {"sample_id", "PCoA1", "PCoA2", "group"}
missing_columns = required_columns.difference(pcoa.columns)

if missing_columns:
    raise ValueError(
        "Missing required columns: " + ", ".join(sorted(missing_columns))
    )

# Colors are assigned explicitly to keep the figure reproducible.
group_colors = {
    "Control": "#377EB8",
    "T2D": "#E41A1C",
}

unknown_groups = sorted(set(pcoa["group"]) - set(group_colors))
if unknown_groups:
    raise ValueError(
        "No color has been defined for: " + ", ".join(unknown_groups)
    )

fig, ax = plt.subplots(figsize=(7, 5.5))

for group, group_data in pcoa.groupby("group", sort=True):
    ax.scatter(
        group_data["PCoA1"],
        group_data["PCoA2"],
        s=110,
        color=group_colors[group],
        edgecolor="black",
        linewidth=0.8,
        label=group,
        zorder=3,
    )

for _, row in pcoa.iterrows():
    ax.annotate(
        row["sample_id"],
        (row["PCoA1"], row["PCoA2"]),
        xytext=(6, 6),
        textcoords="offset points",
        fontsize=10,
    )

ax.axhline(0, color="grey", linewidth=0.8, linestyle="--", zorder=1)
ax.axvline(0, color="grey", linewidth=0.8, linestyle="--", zorder=1)
ax.set_xlabel("PCoA1 (88.92%)")
ax.set_ylabel("PCoA2 (11.08%)")
ax.set_title("PCoA based on Bray–Curtis dissimilarity")
ax.legend(title="Group", frameon=False)
ax.grid(alpha=0.15)

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
fig.tight_layout()
fig.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Figure saved to: {OUTPUT_FILE}")
