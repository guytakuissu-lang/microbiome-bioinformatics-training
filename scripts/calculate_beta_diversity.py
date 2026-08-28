import numpy as np
import pandas as pd


# Load data
feature_table = pd.read_csv(
    "data/example_feature_table.tsv",
    sep="\t",
)

metadata = pd.read_csv(
    "data/example_metadata.tsv",
    sep="\t",
)

# Transform table: samples become rows and ASVs become columns
sample_table = feature_table.set_index("Feature_ID").transpose()
sample_table.index.name = "sample_id"

# Verify that no sample has a total count of zero
library_sizes = sample_table.sum(axis=1)

if (library_sizes == 0).any():
    zero_samples = library_sizes[library_sizes == 0].index.tolist()
    raise ValueError(f"Samples with zero total counts: {zero_samples}")

# Convert counts to proportions
relative_abundance = sample_table.div(library_sizes, axis=0)

sample_ids = relative_abundance.index.tolist()
number_of_samples = len(sample_ids)

# Calculate Bray-Curtis dissimilarities
distance_matrix = np.zeros((number_of_samples, number_of_samples))

for i in range(number_of_samples):
    for j in range(number_of_samples):
        sample_i = relative_abundance.iloc[i].to_numpy()
        sample_j = relative_abundance.iloc[j].to_numpy()

        numerator = np.abs(sample_i - sample_j).sum()
        denominator = (sample_i + sample_j).sum()

        distance_matrix[i, j] = numerator / denominator

bray_curtis = pd.DataFrame(
    distance_matrix,
    index=sample_ids,
    columns=sample_ids,
)
bray_curtis.index.name = "sample_id"

# Verify the distance matrix
if not np.allclose(distance_matrix, distance_matrix.T):
    raise ValueError("The distance matrix is not symmetric.")

if not np.allclose(np.diag(distance_matrix), 0):
    raise ValueError("The matrix diagonal is not zero.")

if np.any(distance_matrix < 0) or np.any(distance_matrix > 1):
    raise ValueError("Bray-Curtis values must be between 0 and 1.")

# Principal Coordinates Analysis using classical scaling
n = number_of_samples
identity = np.eye(n)
centering_matrix = identity - np.ones((n, n)) / n

double_centered = (
    -0.5
    * centering_matrix
    @ (distance_matrix**2)
    @ centering_matrix
)

eigenvalues, eigenvectors = np.linalg.eigh(double_centered)

# Sort eigenvalues and eigenvectors in descending order
order = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[order]
eigenvectors = eigenvectors[:, order]

# Keep only positive eigenvalues
positive = eigenvalues > 1e-12
positive_eigenvalues = eigenvalues[positive]
positive_eigenvectors = eigenvectors[:, positive]

coordinates = positive_eigenvectors * np.sqrt(positive_eigenvalues)

explained_variance = (
    positive_eigenvalues / positive_eigenvalues.sum() * 100
)

pcoa = pd.DataFrame(
    coordinates[:, :2],
    index=sample_ids,
    columns=["PCoA1", "PCoA2"],
)

pcoa.index.name = "sample_id"
pcoa = pcoa.reset_index()

# Add metadata
pcoa = pcoa.merge(
    metadata,
    on="sample_id",
    how="left",
    validate="one_to_one",
)

# Print results
print("======================================")
print("BRAY-CURTIS DISTANCE MATRIX")
print("======================================")
print()
print(bray_curtis.round(3).to_string())

print()
print("======================================")
print("PCoA COORDINATES")
print("======================================")
print()
print(pcoa.round(4).to_string(index=False))

print()
print("Explained variance:")
print(f"PCoA1: {explained_variance[0]:.2f}%")
print(f"PCoA2: {explained_variance[1]:.2f}%")

# Save outputs
bray_curtis.to_csv(
    "results/bray_curtis_distance.tsv",
    sep="\t",
)

pcoa.to_csv(
    "results/pcoa_coordinates.tsv",
    sep="\t",
    index=False,
)

print()
print("Results saved to:")
print("results/bray_curtis_distance.tsv")
print("results/pcoa_coordinates.tsv")
