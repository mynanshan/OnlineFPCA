#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=ci-simu2d
#SBATCH --output=logs/ci-simu2d_%j.out
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument(
  "--seed",
  type = "integer",
  default = 1234,
  help = "Seed used for this experiment."
)
args <- parser$parse_args()
seed <- as.integer(args[["seed"]])

library(fda)
library(Matrix)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/onlineFDAlocalpoly.R")
source("./R/fpcaReg.R")
source("./data_generation/generator.R")

mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

noise_sd <- 0.1
nBatch <- 5
nParams <- 6
nPass <- 5
nBlock <- 100
Ninit <- 100
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize0 <- 0.1
stepsize.min <- 1e-4
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune

exprmt <- "ci2d"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) {
  dir.create(dirpath, recursive = TRUE)
}

n_digit_seed = ceiling(log10(seed + 1))
seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")

set.seed(seed)