from pathlib import Path

import pandas as pd


OUTPUT_PATH = Path("results/compositionality_demo.tsv")


# Toy example: only Taxon_A doubles in absolute abundance.
# Taxon_B and Taxon_C remain biologically unchanged in this scenario.
counts = pd.DataFrame(
    {
        "baseline_count": [100, 100, 100],
        "after_A_doubles_count": [200, 100, 100],
    },
    index=["Taxon_A", "Taxon_B", "Taxon_C"],
)

counts.index.name = "taxon"

# Convert each sample/condition to a closed composition summing to 100%.
relative_abundance = counts.div(counts.sum(axis=0), axis=1) * 100

results = counts.copy()
results["baseline_relative_pct"] = relative_abundance["baseline_count"]
results["after_A_doubles_relative_pct"] = relative_abundance[
    "after_A_doubles_count"
]
results["absolute_count_change"] = (
    results["after_A_doubles_count"] - results["baseline_count"]
)
results["relative_percentage_point_change"] = (
    results["after_A_doubles_relative_pct"]
    - results["baseline_relative_pct"]
)

# Internal checks make the intended lesson explicit and reproducible.
if results.loc["Taxon_B", "absolute_count_change"] != 0:
    raise ValueError("Taxon_B should be unchanged in absolute abundance.")

if results.loc["Taxon_C", "absolute_count_change"] != 0:
    raise ValueError("Taxon_C should be unchanged in absolute abundance.")

if results.loc["Taxon_B", "relative_percentage_point_change"] >= 0:
    raise ValueError("Taxon_B should appear to decrease relatively.")

if results.loc["Taxon_C", "relative_percentage_point_change"] >= 0:
    raise ValueError("Taxon_C should appear to decrease relatively.")

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(OUTPUT_PATH, sep="\t", float_format="%.2f")

print("======================================")
print("COMPOSITIONALITY DEMONSTRATION")
print("======================================")
print()
print(results.round(2).to_string())

print()
print("Relative-abundance column sums:")
print(
    relative_abundance.sum(axis=0).round(2).to_string()
)

print()
print("Interpretation:")
print("- Taxon_A doubles from 100 to 200 counts.")
print("- Taxon_B and Taxon_C remain at 100 counts.")
print("- Nevertheless, B and C fall from 33.33% to 25.00%.")
print("- Their apparent relative decrease is caused by closure to 100%.")

print()
print(f"Results saved to: {OUTPUT_PATH}")
