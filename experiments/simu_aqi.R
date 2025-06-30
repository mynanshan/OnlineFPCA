#!/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/r/4.5.0/bin/Rscript
#SBATCH --job-name=simu_aqi        # Job name
#SBATCH --output=logs/simu_gfr_%j.out         # Standard output file (%j expands to jobID)
#SBATCH --time=3:00:00                   # Maximum runtime (hh:mm:ss)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G                          # Memory allocation

parser <- argparse::ArgumentParser(
  description = "Determine the settings for this simple simulation."
)
parser$add_argument("--seed",
  type = "integer", default = 1234,
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

# load mOpCov codes
mOpCov_path <- "external_codes/mOpCov/"
source(paste0(mOpCov_path, "mOpCov_prep.R"))
Rcpp::sourceCpp(paste0(mOpCov_path, "mOpCov_cpp.cpp"))

# settings
noise <- 0.1
nBatch <- 5
nParams <- 6
nPass <- 5
nBlock <- 100
Ninit <- 100
N <- 5000
nIters.1pass <- seq(nBlock,N,nBlock)
nIters <- seq(nBlock,nPass*N,nBlock)
nRecord.1pass <- round(N/nBlock)
nRecord <- length(nIters)
stepsize <- 2e-1
stepsize.min <- 1e-2
sgd.step.scale <- 1 # when < 1, use a smaller step size for sgd 
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

exprmt <- "aqi"
dirpath <- file.path("experiments", exprmt)
if (!dir.exists(dirpath)) dir.create(dirpath, recursive = TRUE)

n_digit_seed = ceiling(log10(seed))
seedtext <- paste0(rep("0", 4 - n_digit_seed), seed)
filename <- paste0("simu_gfr_sd",seedtext,".csv")