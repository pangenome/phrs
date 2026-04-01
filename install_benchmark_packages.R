#!/usr/bin/env Rscript

# Install required packages for phyper benchmarking
cat("Installing required R packages for benchmarking...\n")

required_packages <- c("microbenchmark", "pryr", "ggplot2", "dplyr", "tibble")

# Install packages that aren't already installed
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
  } else {
    cat("Package", pkg, "already installed\n")
  }
}

cat("All required packages installed successfully!\n")