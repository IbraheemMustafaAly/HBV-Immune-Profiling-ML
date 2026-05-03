# ============================================================
#  R Package Requirements
#  HBV Immune Profiling Study
#  Run this script once to install all required packages
# ============================================================

# Check R version
cat("R version required: >= 4.3.0\n")
cat("Current R version:", R.version.string, "\n\n")

# Core packages
required_packages <- c(
  # Data manipulation
  "tidyverse",    # 2.0.0  — ggplot2, dplyr, tidyr, readr, purrr
  "dplyr",        # 1.1.x
  "tidyr",        # 1.3.x
  "readr",        # 2.1.x

  # Statistical analysis
  "rstatix",      # 0.7.2  — pipe-friendly stats
  "HardyWeinberg",# 1.7.5  — HWE equilibrium testing
  "epitools",     # 0.5-10 — odds ratios

  # Visualization
  "ggplot2",      # 3.4.x  — base plotting (via tidyverse)
  "ggpubr",       # 0.6.0  — publication-ready figures with significance
  "ggrepel",      # 0.9.x  — text label repulsion

  # Multivariate analysis
  "FactoMineR",   # 2.9    — PCA
  "factoextra",   # 1.0.7  — PCA visualization

  # Correlation & heatmaps
  "corrplot",     # 0.92   — correlation heatmaps
  "pheatmap",     # 1.0.12 — hierarchical clustering heatmaps
  "ComplexHeatmap", # 2.x  — advanced heatmaps (Bioconductor)
  "circlize",     # 0.4.x  — color functions for heatmaps

  # Color palettes
  "RColorBrewer", # 1.1-3  — color palettes
  "scales",       # 1.2.x  — scale functions

  # Figure composition
  "patchwork"     # 1.2.0  — multi-panel figure assembly
)

# Install CRAN packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing:", pkg))
    install.packages(pkg, dependencies = TRUE,
                     repos = "https://cloud.r-project.org")
  } else {
    message(paste("Already installed:", pkg,
                  "-", packageVersion(pkg)))
  }
}

lapply(required_packages[required_packages != "ComplexHeatmap"],
       install_if_missing)

# Install ComplexHeatmap from Bioconductor
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}

# Print session info for reproducibility
cat("\n============================================\n")
cat("Package versions installed:\n")
cat("============================================\n")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-20s %s\n", pkg,
                as.character(packageVersion(pkg))))
  }
}

cat("\nR Session Info:\n")
print(sessionInfo())
