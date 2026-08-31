suppressPackageStartupMessages({
  library(ANCOMBC)
  library(BiocManager)
})

cat("======================================\n")
cat("ANCOM-BC2 ENVIRONMENT CHECK\n")
cat("======================================\n\n")

cat("R version:", R.version.string, "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
cat("ANCOMBC version:", as.character(packageVersion("ANCOMBC")), "\n")
cat("ANCOMBC location:", find.package("ANCOMBC"), "\n")

ancombc2_available <- exists(
  "ancombc2",
  where = asNamespace("ANCOMBC"),
  inherits = FALSE
)

cat("ancombc2 function available:", ancombc2_available, "\n")

if (!ancombc2_available) {
  stop("The installed ANCOMBC package does not expose ancombc2().")
}

cat("\nProject library paths:\n")
print(.libPaths())

cat("\nSession information:\n")
print(sessionInfo())
