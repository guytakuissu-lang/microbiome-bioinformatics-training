import pandas as pd

metadata = pd.read_csv("data/example_metadata.tsv", sep="\t")

print("Metadata categorical validation")
print("===============================")

# Allowed values
allowed_groups = {"T2D", "Control"}
allowed_sex = {"M", "F"}

# Check group
invalid_groups = set(metadata["group"]) - allowed_groups

print("\nGroup validation:")
if invalid_groups:
    print("ERROR: Invalid group values:", invalid_groups)
else:
    print("PASS: All group values are valid.")

# Check sex
invalid_sex = set(metadata["sex"]) - allowed_sex

print("\nSex validation:")
if invalid_sex:
    print("ERROR: Invalid sex values:", invalid_sex)
else:
    print("PASS: All sex values are valid.")

# Check age
invalid_age = metadata[
    (metadata["age"] < 18) | (metadata["age"] > 120)
]

print("\nAge validation:")
if len(invalid_age) > 0:
    print("ERROR: Invalid age values detected.")
    print(invalid_age)
else:
    print("PASS: All age values are within the expected range.")
