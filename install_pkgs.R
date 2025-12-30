install.packages("pak")
library(pak)
required_pkgs <- c(
  "fda", "face", "fdapace",
  "Metrics", "fastmatrix", "RSpectra",
  "mgcv", "argparse", "Rcpp", "rTensor", "gslnls",
  "stringr", "tidyr", "dplyr", "readr", "cowplot",
  "sf", "spData", "giscoR",
  "lubridate", "readxl", "Rdimtools", "sm", 
  "foreach", "doFuture", "ManifoldOptim", "remotes"
)
# NOTE: date 2025-06-29
# At present time, the newest RcppArnadillo is unstable
# and it could make mOpCov fail, giving "index out of bound" error.
# Restore to an order version:
pkg_install(required_pkgs)
remotes::install_version("RcppArmadillo", version = "0.12.8.4.0")
install.packages(
  "external_codes/fpca_0.2-1.tar.gz", repos = NULL,
  type="source", dependencies = TRUE
)