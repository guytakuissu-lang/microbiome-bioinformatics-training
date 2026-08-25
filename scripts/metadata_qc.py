import pandas as pd
import sys

METADATA_FILE = "data/example_metadata.tsv"

print("======================================")
print("MICROBIOME METADATA QC")
print("======================================")
print(f"File: {METADATA_FILE}")

# Load metadata
try:
    metadata = pd.read_csv(METADATA_FILE, sep="\t")
except Exception as e:
    print(f"ERROR: Could not read metadata file: {e}")
    sys.exit(1)

errors = 0

# -----------------------------
# Basic structure
# -----------------------------

print("\nBasic structure")
print("----------------")

rows, columns = metadata.shape

print(f"Samples: {rows}")
print(f"Variables: {columns}")

required_columns = {"sample_id", "group", "age", "sex"}

missing_columns = required_columns - set(metadata.columns)

if missing_columns:
    print(f"ERROR: Missing required columns: {missing_columns}")
    errors += 1
else:
    print("PASS: Required columns are present.")

# -----------------------------
# Missing values
# -----------------------------

print("\nMissing values")
print("----------------")

missing = metadata.isna().sum()

if missing.sum() == 0:
    print("PASS: No missing values detected.")
else:
    print("ERROR: Missing values detected.")
    print(missing[missing > 0])
    errors += 1

# -----------------------------
# Sample IDs
# -----------------------------

print("\nSample ID validation")
print("---------------------")

if metadata["sample_id"].is_unique:
    print("PASS: All sample IDs are unique.")
else:
    print("ERROR: Duplicate sample IDs detected.")
    errors += 1

empty_ids = (metadata["sample_id"].astype(str).str.strip() == "").sum()

if empty_ids == 0:
    print("PASS: No empty sample IDs.")
else:
    print(f"ERROR: {empty_ids} empty sample IDs detected.")
    errors += 1

# -----------------------------
# Group validation
# -----------------------------

print("\nGroup validation")
print("----------------")

allowed_groups = {"T2D", "Control"}
invalid_groups = set(metadata["group"]) - allowed_groups

if not invalid_groups:
    print("PASS: All group values are valid.")
else:
    print(f"ERROR: Invalid group values: {invalid_groups}")
    errors += 1

# -----------------------------
# Sex validation
# -----------------------------

print("\nSex validation")
print("----------------")

allowed_sex = {"M", "F"}
invalid_sex = set(metadata["sex"]) - allowed_sex

if not invalid_sex:
    print("PASS: All sex values are valid.")
else:
    print(f"ERROR: Invalid sex values: {invalid_sex}")
    errors += 1

# -----------------------------
# Age validation
# -----------------------------

print("\nAge validation")
print("----------------")

invalid_age = metadata[
    (metadata["age"] < 18) | (metadata["age"] > 120)
]

if len(invalid_age) == 0:
    print("PASS: All ages are within the expected range.")
else:
    print("ERROR: Invalid age values detected.")
    print(invalid_age)
    errors += 1

# -----------------------------
# Group distribution
# -----------------------------

print("\nGroup distribution")
print("-------------------")
print(metadata["group"].value_counts())

# -----------------------------
# Sex distribution
# -----------------------------

print("\nSex distribution")
print("----------------")
print(metadata["sex"].value_counts())

# -----------------------------
# Age summary
# -----------------------------

print("\nAge summary")
print("-----------")
print(metadata["age"].describe())

# -----------------------------
# Final status
# -----------------------------

print("\n======================================")

if errors == 0:
    print("STATUS: PASS")
else:
    print(f"STATUS: FAIL ({errors} error(s))")

print("======================================")

sys.exit(0 if errors == 0 else 1)
