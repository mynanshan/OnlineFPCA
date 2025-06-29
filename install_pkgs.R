install.packages("pak")
library(pak)
required_pkgs <- c(
  "fda", "face", "fdapace",
  "Metrics", "fastmatrix", "RSpectra",
  "mgcv", "argparse", "Rcpp",
  "stringr", "tidyr", "dplyr", "readr",
  "lubridate", "readxl", "Rdimtools", "sm"
)
pkg_install(required_pkgs)
install.packages(
  "external_codes/fpca_0.2-1.tar.gz", repos = NULL,
  type="source", dependencies = TRUE
)