from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


METADATA_PATH = Path("data/example_metadata.tsv")
RESULTS_DIR = Path("results")
OUTPUT_PATH = Path("figures/taxonomic_composition.png")

PHYLUM_COLORS = {
    "Actinomycetota": "#E69F00",
    "Bacillota": "#56B4E9",
    "Bacteroidota": "#009E73",
    "Pseudomonadota": "#D55E00",
    "Verrucomicrobiota": "#CC79A7",
}

GENUS_COLORS = {
    "Akkermansia": "#8DD3C7",
    "Bacteroides": "#FFFFB3",
    "Bifidobacterium": "#BEBADA",
    "Blautia": "#FB8072",
    "Escherichia-Shigella": "#80B1D3",
    "Faecalibacterium": "#FDB462",
}


metadata = pd.read_csv(METADATA_PATH, sep="\t")

if not {"sample_id", "group"}.issubset(metadata.columns):
    raise ValueError("Metadata must contain sample_id and group columns.")

if metadata["sample_id"].duplicated().any():
    raise ValueError("Metadata contain duplicated sample_id values.")

sample_order = metadata["sample_id"].tolist()
group_by_sample = metadata.set_index("sample_id")["group"]
sample_labels = [
    f"{sample_id}\n{group_by_sample[sample_id]}"
    for sample_id in sample_order
]


def load_relative_abundance(rank: str) -> pd.DataFrame:
    """Load and validate a taxonomic relative-abundance table."""
    path = RESULTS_DIR / f"{rank.lower()}_relative_abundance.tsv"
    table = pd.read_csv(path, sep="\t", index_col=0)

    missing_samples = sorted(set(sample_order) - set(table.columns))
    extra_samples = sorted(set(table.columns) - set(sample_order))

    if missing_samples or extra_samples:
        raise ValueError(
            f"{rank} sample mismatch. Missing: {missing_samples}; "
            f"extra: {extra_samples}"
        )

    table = table[sample_order]

    if table.isna().any().any():
        raise ValueError(f"{rank} table contains missing values.")

    if (table < 0).any().any():
        raise ValueError(f"{rank} table contains negative values.")

    if not table.sum(axis=0).round(6).eq(100).all():
        raise ValueError(f"{rank} relative abundances do not sum to 100%.")

    return table


def plot_stacked_composition(
    ax: plt.Axes,
    table: pd.DataFrame,
    colors: dict[str, str],
    title: str,
) -> None:
    """Draw one stacked taxonomic-composition panel."""
    undefined_colors = sorted(set(table.index) - set(colors))
    if undefined_colors:
        raise ValueError(
            "No color defined for: " + ", ".join(undefined_colors)
        )

    bottom = pd.Series(0.0, index=sample_order)

    for taxon in table.index:
        values = table.loc[taxon]
        ax.bar(
            sample_labels,
            values,
            bottom=bottom,
            label=taxon,
            color=colors[taxon],
            edgecolor="white",
            linewidth=0.5,
            width=0.72,
        )
        bottom += values

    ax.set_title(title)
    ax.set_ylabel("Relative abundance (%)")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", alpha=0.2)
    ax.set_axisbelow(True)
    ax.legend(
        title="Taxon",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )


phylum = load_relative_abundance("Phylum")
genus = load_relative_abundance("Genus")

fig, axes = plt.subplots(2, 1, figsize=(10, 11), sharex=False)

plot_stacked_composition(
    axes[0],
    phylum,
    PHYLUM_COLORS,
    "A. Composition at phylum level",
)

plot_stacked_composition(
    axes[1],
    genus,
    GENUS_COLORS,
    "B. Composition at genus level",
)

fig.suptitle(
    "Taxonomic composition of simulated microbiome samples",
    fontsize=15,
    y=0.995,
)
fig.text(
    0.5,
    0.005,
    "Simulated data for bioinformatics training; no biological inference.",
    ha="center",
    fontsize=9,
    color="dimgray",
)

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
fig.tight_layout(rect=(0, 0.025, 1, 0.98))
fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Figure saved to: {OUTPUT_PATH}")
