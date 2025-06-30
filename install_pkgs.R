install.packages("pak")
library(pak)
required_pkgs <- c(
  "fda", "face", "fdapace",
  "Metrics", "fastmatrix", "RSpectra",
  "mgcv", "argparse", "Rcpp",
  "stringr", "tidyr", "dplyr", "readr",
  "lubridate", "readxl", "Rdimtools", "sm"
)
# NOTE: date 2025-06-29
# At present time, the newest RcppArnadillo is unstable
# and it could make mOpCov fail, giving index out of bound error.
# Restore to an order version:
# remotes::install_version("RcppArmadillo", version = "0.12.8.4.0")
pkg_install(required_pkgs)
install.packages(
  "external_codes/fpca_0.2-1.tar.gz", repos = NULL,
  type="source", dependencies = TRUE
)