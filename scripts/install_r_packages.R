cran_repo <- "https://cloud.r-project.org"

ensure_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = cran_repo)
  }
}

cran_packages <- c(
  "Seurat",
  "patchwork",
  "ggplot2",
  "dplyr",
  "readr",
  "tibble",
  "jsonlite",
  "hdf5r"
)

for (pkg in cran_packages) {
  ensure_package(pkg)
}
