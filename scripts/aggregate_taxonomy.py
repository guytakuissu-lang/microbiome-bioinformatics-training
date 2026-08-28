from pathlib import Path

import pandas as pd


FEATURE_TABLE_PATH = Path("data/example_feature_table.tsv")
TAXONOMY_PATH = Path("data/example_taxonomy.tsv")
RESULTS_DIR = Path("results")


feature_table = pd.read_csv(FEATURE_TABLE_PATH, sep="\t")
taxonomy = pd.read_csv(TAXONOMY_PATH, sep="\t", dtype=str)

if "Feature_ID" not in feature_table.columns:
    raise ValueError("Feature table does not contain Feature_ID.")

if "Feature_ID" not in taxonomy.columns:
    raise ValueError("Taxonomy table does not contain Feature_ID.")

sample_ids = feature_table.columns.drop("Feature_ID").tolist()

# Confirm that abundance values are numeric and biologically valid counts.
abundance = feature_table[sample_ids].apply(pd.to_numeric, errors="raise")

if abundance.isna().any().any():
    raise ValueError("The feature table contains missing abundance values.")

if (abundance < 0).any().any():
    raise ValueError("The feature table contains negative abundance values.")

if (abundance.sum(axis=0) == 0).any():
    zero_samples = abundance.columns[abundance.sum(axis=0) == 0].tolist()
    raise ValueError(f"Samples with zero total counts: {zero_samples}")

# Merge each ASV with its taxonomy. One-to-one validation prevents silent
# duplication of counts when taxonomy identifiers are duplicated.
merged = feature_table.merge(
    taxonomy,
    on="Feature_ID",
    how="left",
    validate="one_to_one",
)


def aggregate_by_rank(rank: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate raw counts and relative abundance at one taxonomic rank."""
    if rank not in merged.columns:
        raise ValueError(f"Taxonomic rank not found: {rank}")

    if merged[rank].isna().any():
        missing_features = merged.loc[
            merged[rank].isna(), "Feature_ID"
        ].tolist()
        raise ValueError(
            f"Missing {rank} assignment for: {missing_features}"
        )

    counts = merged.groupby(rank, sort=True)[sample_ids].sum()
    relative_abundance = counts.div(counts.sum(axis=0), axis=1) * 100

    return counts, relative_abundance


RESULTS_DIR.mkdir(parents=True, exist_ok=True)

for rank in ["Phylum", "Genus"]:
    counts, relative_abundance = aggregate_by_rank(rank)
    rank_name = rank.lower()

    counts_path = RESULTS_DIR / f"{rank_name}_counts.tsv"
    relative_path = RESULTS_DIR / f"{rank_name}_relative_abundance.tsv"

    counts.to_csv(counts_path, sep="\t")
    relative_abundance.to_csv(relative_path, sep="\t")

    print("======================================")
    print(f"{rank.upper()} RELATIVE ABUNDANCE (%)")
    print("======================================")
    print(relative_abundance.round(2).to_string())

    print()
    print("Column sums (%):")
    print(relative_abundance.sum(axis=0).round(2).to_string())

    print()
    print(f"Saved: {counts_path}")
    print(f"Saved: {relative_path}")
    print()

print("Taxonomic aggregation completed successfully.")
