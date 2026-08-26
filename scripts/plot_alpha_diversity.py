import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Load alpha-diversity results
alpha = pd.read_csv(
    "results/alpha_diversity.tsv",
    sep="\t"
)

# Colours assigned to clinical groups
group_colours = {
    "T2D": "#D55E00",
    "Control": "#0072B2"
}

colours = alpha["group"].map(group_colours)

metrics = [
    ("observed_asvs", "Observed ASVs"),
    ("shannon", "Shannon diversity"),
    ("simpson_1_minus_D", "Simpson diversity (1-D)")
]

fig, axes = plt.subplots(
    1,
    3,
    figsize=(12, 4.5)
)

for ax, (column, title) in zip(axes, metrics):

    ax.bar(
        alpha["sample_id"],
        alpha[column],
        color=colours,
        edgecolor="black",
        linewidth=0.8
    )

    ax.set_title(title)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Index value")

    # Add values above bars
    for position, value in enumerate(alpha[column]):
        ax.text(
            position,
            value,
            f"{value:.3f}" if column != "observed_asvs" else f"{int(value)}",
            ha="center",
            va="bottom",
            fontsize=9
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

legend_elements = [
    Patch(
        facecolor=group_colours["T2D"],
        edgecolor="black",
        label="T2D"
    ),
    Patch(
        facecolor=group_colours["Control"],
        edgecolor="black",
        label="Control"
    )
]

fig.legend(
    handles=legend_elements,
    title="Group",
    loc="upper center",
    ncol=2,
    frameon=False
)

fig.suptitle(
    "Alpha diversity of simulated microbiome samples",
    y=1.05
)

plt.tight_layout()

output_path = "results/alpha_diversity.png"

plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Alpha-diversity figure saved to:")
print(output_path)
