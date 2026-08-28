import pandas as pd


TAXONOMY_PATH = "data/example_taxonomy.tsv"
FEATURE_TABLE_PATH = "data/example_feature_table.tsv"

REQUIRED_TAXONOMY_COLUMNS = [
    "Feature_ID",
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
]


# Load all taxonomy fields as text so that their integrity can be checked.
taxonomy = pd.read_csv(TAXONOMY_PATH, sep="\t", dtype=str)
feature_table = pd.read_csv(FEATURE_TABLE_PATH, sep="\t")

print("======================================")
print("TAXONOMY QUALITY CONTROL")
print("======================================")

# Check the expected structure.
missing_columns = [
    column
    for column in REQUIRED_TAXONOMY_COLUMNS
    if column not in taxonomy.columns
]

if missing_columns:
    raise ValueError(
        "Missing taxonomy columns: " + ", ".join(missing_columns)
    )

if "Feature_ID" not in feature_table.columns:
    raise ValueError("Feature table does not contain a Feature_ID column.")

# Restrict subsequent checks to the expected taxonomy columns.
taxonomy = taxonomy[REQUIRED_TAXONOMY_COLUMNS].copy()

# Remove accidental whitespace around identifiers and taxonomic labels.
for column in REQUIRED_TAXONOMY_COLUMNS:
    taxonomy[column] = taxonomy[column].str.strip()

# Convert empty strings to missing values.
taxonomy = taxonomy.replace("", pd.NA)

duplicate_taxonomy_ids = taxonomy.loc[
    taxonomy["Feature_ID"].duplicated(keep=False), "Feature_ID"
].dropna().unique().tolist()

feature_ids = set(feature_table["Feature_ID"].astype(str))
taxonomy_ids = set(taxonomy["Feature_ID"].dropna())

features_without_taxonomy = sorted(feature_ids - taxonomy_ids)
taxonomy_without_features = sorted(taxonomy_ids - feature_ids)

missing_values = taxonomy.isna().sum()

print()
print(f"Taxonomy dimensions: {taxonomy.shape}")
print(f"Number of feature IDs: {taxonomy['Feature_ID'].nunique()}")
print(f"Duplicate feature IDs: {duplicate_taxonomy_ids}")

print()
print("Missing values by column:")
print(missing_values.to_string())

print()
print("Features without taxonomy:")
print(features_without_taxonomy)

print()
print("Taxonomy IDs absent from feature table:")
print(taxonomy_without_features)

print()
print("Number of unique taxa by rank:")
for rank in REQUIRED_TAXONOMY_COLUMNS[1:]:
    print(f"{rank}: {taxonomy[rank].nunique()}")

# Stop the workflow when a critical integrity problem is detected.
problems = []

if duplicate_taxonomy_ids:
    problems.append("duplicate taxonomy Feature_ID values")
if missing_values.sum() > 0:
    problems.append("missing taxonomy values")
if features_without_taxonomy:
    problems.append("features without taxonomy")
if taxonomy_without_features:
    problems.append("taxonomy IDs absent from feature table")

print()
if problems:
    raise ValueError("Taxonomy QC failed: " + "; ".join(problems))

print("All taxonomy quality-control checks passed.")
print("======================================")
