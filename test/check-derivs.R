library(fda)
library(Matrix)
library(foreach)

source("./R/fdaMdim.R")
source("./R/helper.R")
source("./R/manifoldUtils.R")
source("./R/kernelUtils.R")
source("./R/onlineFPCA.R")
source("./R/onlineFDAlocalpoly.R")
source("./R/fpcaReg.R")
source("./data_generation/generator.R")

noise_sd <- 0.5
nBatch <- 5
nParams <- 6
nPass <- 3
nBlock <- 100
Ninit <- 100
initMethod <- "face"
N <- 5000
nIters.1pass <- seq(nBlock, N, nBlock)
nIters <- seq(nBlock, nPass * N, nBlock)
nRecord.1pass <- round(N / nBlock)
nRecord <- length(nIters)
stepsize <- 2e-1
stepsize.min <- 5e-2
sgd.step.scale <- 1 # use a smaller step size for sgd
nRoundNoTune <- 1
nRoundTune <- nRecord.1pass - nRoundNoTune
asgd.use <- TRUE

seed <- 2501
set.seed(seed)

rp <- get_rp.yang2021(alpha = 2, npc = 10)
m_min <- NULL
m_max <- NULL
m_mean <- 6
m_sd <- 2
m_type <- "gaussian"

t0 <- rp$t0
t1 <- rp$t1
D <- rp$ndim

# evaluate true mean and eigenfunctions
neval <- 101
evalGrid <- seq(t0, t1, length.out = neval)
muTrueEval <- rp$meanfun(evalGrid)
PhiTrueEval <- rp$eigfun(evalGrid)
lambdaTrue <- rp$eigval

q <- npc <- 3

dat <- get_measurements(
  rp,
  n = N,
  m_min = m_min,
  m_max = m_max,
  m_mean = m_mean,
  m_sd = m_sd,
  design_type = "random",
  m_type = m_type,
  sigma = noise_sd
)
tgrid <- dat$tgrid

fdata_generator <- function(n, total_count) {
  idx <- (total_count):(total_count + n - 1) %% N + 1
  list(
    Ly = dat$Ly[idx],
    Ltid = dat$Ltid[idx],
    Lt = dat$Lt[idx],
    Lmi = dat$Lmi[idx]
  )
}

nbasis <- 7
basis <- create.bspline.basis(c(t0, t1), nbasis = nbasis, norder = 4)
p <- nbasis

muTrueFunc <- smooth_basis(evalGrid, muTrueEval, basis)
phiTrueFunc <- smooth_basis(evalGrid, PhiTrueEval, basis)
rmseMuBest <- Metrics::rmse(muTrueEval, eval_fd(evalGrid, muTrueFunc))
rmsePhiBest <- Metrics::rmse(PhiTrueEval, eval_fd(evalGrid, phiTrueFunc))

B <- eval_basis(tgrid, basis)
G <- get_basis_inprod_matrix(basis)
Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
sqrtmObj <- pracma::sqrtm(as.matrix(G))
sqrtG <- sqrtmObj$B
sqrtGinv <- sqrtmObj$Binv
invG <- Matrix::solve(G)
GR <- Matrix::chol(G)

tmpdat <- data.frame(
  argvals = unlist(dat$Lt[1:Ninit]),
  subj = rep(1:Ninit, dat$Lmi[1:Ninit]),
  y = unlist(dat$Ly[1:Ninit])
)
fitFace <- face::face.sparse(
  tmpdat,
  center = FALSE,
  argvals.new = evalGrid,
  knots = p
)

muInitEval <- fitFace$mu.new
theta_muInit <- smooth_basis(
  evalGrid,
  muInitEval,
  basis,
  lambda = 1e-10
)$coefs
eig <- RSpectra::eigs_sym(fitFace$Chat.new, q)
PhiInitEval <- flip_direc(eig$vectors, PhiTrueEval[, 1:q])
ThetaInit <- smooth_basis(evalGrid, PhiInitEval, basis, lambda = 1e-10)$coefs
lambdaInit <- eig$values
sigma2Init <- fitFace$sigma2

norm_factor <- diag(t(ThetaInit) %*% G %*% ThetaInit)
ThetaInit <- sweep(ThetaInit, 2, sqrt(norm_factor), "/")
lambdaInit <- lambdaInit * norm_factor
PhiInitEval <- eval_fd(evalGrid, FuncData(ThetaInit, basis))

nBlockIter <- nBlock / nBatch
nIter1pass <- N / nBatch
asgdIterStart <- round(nRecord.1pass * 0.7 * nBlockIter)
adamIterEnd <- round(nRecord.1pass * 1.7 * nBlockIter)

inits <- list(Theta = ThetaInit, lambda = lambdaInit, sigma2 = sigma2Init)

optMethod <- "Adam"

# fit <- fpca.sgd(
#   fdata_generator,
#   tgrid,
#   inits = inits,
#   meanfun = FALSE,
#   tau = NULL,
#   tau.control = list(
#     ntau = nParams,
#     nselect = 2,
#     maxtau = 1e-1,
#     mintau = 1e-6
#   ),
#   nbatch = nBatch,
#   maxIter = nPass * nIter1pass,
#   stepsize = stepsize,
#   stepsize.decayrate = 0.6,
#   stepsize.min = stepsize.min,
#   nIter.slowerdecay = floor(0.8 * nIter1pass),
#   stepsize.decayrate.slow = 0.3,
#   nIter.adam = ifelse(optMethod == "SGD", 0, adamIterEnd),
#   asgd.use = TRUE,
#   asgd.start = asgdIterStart,
#   coord.scaling = FALSE,
#   nIter.1stTune = nRoundNoTune * nBlockIter,
#   nIter.lastTune = nIter1pass,
#   nIter.tauNoIncrease = floor(0.7 * nIter1pass),
#   period.tune = nBlockIter,
#   period.record = nBlockIter,
#   verbose = FALSE
# )


# Check Gradients ------------------------------------------------------

library(nloptr)
N0 <- 20

ThetaTrue = phiTrueFunc$coefs[,1:q]
lambdaTrue
sigma2true = noise_sd^2

Theta = ThetaInit |> manifold.Stiefel.retract(G)
lambda = lambdaInit
sigma2 = sigma2Init

# Theta = ThetaInit
# lambda = lambdaInit
# sigma2 = sigma2Init


compute_grads <- function(
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(length(Ltid) == n)
  grad_Theta <- matrix(0, p, q)
  grad_eta <- numeric(q)
  grad_zeta <- 0

  for (i in seq_len(n)) {
    mi <- length(Ltid[[i]])
    Bi <- B[Ltid[[i]], , drop = F]
    Yi <- Ly[[i]]
    if (!is.null(theta_mu)) {
      Yi <- Yi - Bi %*% theta_mu
    }
    Phi <- as.matrix(Bi %*% Theta)
    PhiT_Yi <- crossprod(Phi, Yi)
    Q <- crossprod(Phi, Phi) + sigma2 * diag(1 / lambda, q, q)
    invQ <- solve(Q)
    Phit_Phi <- crossprod(Phi)
    invQ_PhiT_Yi <- invQ %*% PhiT_Yi

    grad_Theta <- grad_Theta +
      2 / n *
        as.matrix(
          crossprod(Bi, Phi) %*% (
            1 / sigma2 * invQ_PhiT_Yi %*% t(invQ_PhiT_Yi)
              + invQ
          )
            - 1 / sigma2 * crossprod(Bi, Yi) %*% t(invQ_PhiT_Yi)
        )

    grad_eta <- grad_eta +
      1 / n * (
        -invQ_PhiT_Yi^2 / lambda + rep(1, q)
          - sigma2 * diag(invQ) / lambda
      )

    grad_zeta <- grad_zeta +
      1 / n * (
        -sum(Yi^2) / sigma2
          + sum(invQ_PhiT_Yi * PhiT_Yi) / sigma2
          + sum(invQ_PhiT_Yi^2 / lambda)
          + mi - q
          + sigma2 * sum(diag(invQ) / lambda)
      )
  }
  grad_Theta <- grad_Theta + 2 * tau * as.matrix(Omega %*% Theta)
  list(
    grad_Theta = as.matrix(Matrix::solve(G, grad_Theta)),
    grad_eta = as.numeric(grad_eta),
    grad_zeta = as.numeric(grad_zeta)
  )
}


check_grad_Theta <- check.derivatives(
  inits$Theta,
  \(x) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], matrix(x, p, q), inits$lambda, inits$sigma2,
      stats = "loss"
    )[["fval"]]
  },
  \(x) {
    as.vector(
      G %*% compute_grads(
        dat$Ly[1:N0], dat$Ltid[1:N0], matrix(x, p, q), inits$lambda, inits$sigma2
      )[["grad_Theta"]]
    )
  }
)

check_grad_eta <- check.derivatives(
  log(inits$lambda),
  \(x) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], inits$Theta, exp(x), inits$sigma2,
      stats = "loss"
    )[["fval"]]
  },
  \(x) {
    compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], inits$Theta, exp(x), inits$sigma2
    )[["grad_eta"]]
  }
)

check_grad_zeta <- check.derivatives(
  log(inits$sigma2),
  \(x) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], inits$Theta, inits$lambda, exp(x),
      stats = "loss"
    )[["fval"]]
  },
  \(x) {
    compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], inits$Theta, inits$lambda, exp(x)
    )[["grad_zeta"]]
  }
)


# Check Hessian Theta ----------------------------------------------------

Hess_x_Theta_apply <- function(
  Delta,
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(ncol(Delta) == q)
  stopifnot(nrow(Delta) == p)
  stopifnot(length(Ltid) == n)
  H_Theta_Theta <- matrix(0, p, q)
  H_eta_Theta <- numeric(q)
  H_zeta_Theta <- 0

  for (i in seq_len(n)) {
    mi <- length(Ltid[[i]])
    Bi <- B[Ltid[[i]], , drop = F]
    Yi <- Ly[[i]]
    if (!is.null(theta_mu)) {
      Yi <- Yi - Bi %*% theta_mu
    }
    Phi <- as.matrix(Bi %*% Theta)
    Bt_Yi <- as.matrix(crossprod(Bi, Yi))
    PhiT_Yi <- crossprod(Phi, Yi)
    Phit_Phi <- crossprod(Phi)
    Q <- crossprod(Phi, Phi) + sigma2 * diag(1 / lambda, q, q)
    invQ <- solve(Q)
    Theta_invQ <- Theta %*% invQ
    invQ_PhiT_Yi <- invQ %*% PhiT_Yi
    BDel <- as.matrix(Bi %*% Delta)
    BDelT_Yi <- crossprod(BDel, Yi)
    Phit_BDel <- crossprod(Phi, BDel)
    BtB <- as.matrix(crossprod(Bi))
    Delta_invQ <- Delta %*% invQ
    invQ_BDelT_Yi <- invQ %*% BDelT_Yi

    ### H_Theta_Theta

    term1 <- 2 / sigma2 * BtB %*% (
      Delta_invQ %*% t(Theta) + Theta_invQ %*% t(Delta)
      - Theta_invQ %*% (Phit_BDel + t(Phit_BDel)) %*% t(Theta_invQ) 
    ) %*% Bt_Yi %*% t(invQ_PhiT_Yi)
    term2 <- 2 / sigma2 * (BtB %*% Theta_invQ %*% t(Theta) - diag(1,p,p)) %*%
      Bt_Yi %*% (t(Bt_Yi) %*% (
        Delta_invQ - Theta_invQ %*% (Phit_BDel + t(Phit_BDel)) %*% invQ
      ))
    term3 <- 2 * BtB %*% {
      Delta - Theta_invQ %*% (Phit_BDel + t(Phit_BDel))
    } %*% invQ
    term4 <- 2 * tau * as.matrix(Omega) %*% Delta

    H_Theta_Theta <- H_Theta_Theta +
      1 / n * (term1 + term2 + term3 + term4)

    ### H_eta_Theta

    term1 <- -2 * (
      invQ_PhiT_Yi * (
        - invQ %*% (Phit_BDel + t(Phit_BDel)) %*% invQ_PhiT_Yi
        + invQ_BDelT_Yi
      ) / lambda
    ) |> as.vector()
    term2 <- sigma2 / lambda * diag(
      invQ %*% (Phit_BDel + t(Phit_BDel)) %*% invQ
    )

    H_eta_Theta <- H_eta_Theta + 1 / n * (term1 + term2)

    ### H_zeta_Theta

    term1 <- 2 / sigma2 * sum(BDelT_Yi * invQ_PhiT_Yi)
    term2 <- -1 / sigma2 * sum(
      invQ_PhiT_Yi * (Phit_BDel + t(Phit_BDel)) %*% invQ_PhiT_Yi
    )
    term3 <- 2 * sum(invQ_BDelT_Yi * invQ_PhiT_Yi / lambda)
    term4 <- -2 * sum(
      invQ_PhiT_Yi * (Phit_BDel + t(Phit_BDel)) %*%
        invQ %*% (invQ_PhiT_Yi / lambda)
    )
    term5 <- -sigma2 * sum(
      invQ * (Phit_BDel + t(Phit_BDel)) %*% invQ %*% diag(1/lambda, q, q)
    )

    H_zeta_Theta <- H_zeta_Theta + 1 / n * (
      term1 + term2 + term3 + term4 + term5
    )
  }

  list(
    H_Theta_Theta = as.matrix(invG %*% H_Theta_Theta),
    H_eta_Theta = as.vector(H_eta_Theta),
    H_zeta_Theta = as.numeric(H_zeta_Theta)
  )
}

Delta <- matrix(rnorm(p*q), p, q)
Delta <- manifold.Stiefel.project(Delta, Theta, G)
Hess_x_Theta_apply(Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2)

check_H_Theta_Theta <- check.derivatives(
  0,
  \(h) {
    compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta + h * Delta, lambda, sigma2
    )[["grad_Theta"]] |> as.vector()
  },
  \(h) {
    Hess_x_Theta_apply(
      Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_Theta_Theta"]] |> as.vector()
  }
)

check_H_eta_Theta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta + h * Delta, lambda, sigma2
    )[["grad_eta"]]
  },
  \(h) {
    Hess_x_Theta_apply(
      Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_eta_Theta"]]
  }
)

check_H_zeta_Theta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta + h * Delta, lambda, sigma2
    )[["grad_zeta"]]
  },
  \(h) {
    Hess_x_Theta_apply(
      Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_zeta_Theta"]]
  }
)


# Check Hessian eta ------------------------------------------------------

Hess_x_eta_apply <- function(
  delta,
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(length(delta) == q)
  stopifnot(length(Ltid) == n)
  H_Theta_eta <- matrix(0, p, q)
  H_eta_eta <- numeric(q)
  H_zeta_eta <- 0

  for (i in seq_len(n)) {
    mi <- length(Ltid[[i]])
    Bi <- B[Ltid[[i]], , drop = F]
    Yi <- Ly[[i]]
    if (!is.null(theta_mu)) {
      Yi <- Yi - Bi %*% theta_mu
    }
    BtB <- as.matrix(crossprod(Bi))
    Phi <- as.matrix(Bi %*% Theta)
    Bt_Yi <- as.matrix(crossprod(Bi, Yi))
    PhiT_Yi <- crossprod(Phi, Yi)
    Phit_Phi <- crossprod(Phi)
    Q <- crossprod(Phi, Phi) + sigma2 * diag(1 / lambda, q, q)
    invQ <- solve(Q)
    Theta_invQ <- Theta %*% invQ
    invQ_PhiT_Yi <- invQ %*% PhiT_Yi

    ### H_Theta_eta

    term1 <- 2 * (
      BtB %*% Theta_invQ %*% diag(delta/lambda, q, q) %*%
        invQ_PhiT_Yi %*% t(invQ_PhiT_Yi)
    )
    term2 <- 2 * (BtB %*% Theta_invQ %*% t(Theta) - diag(1,p,p)) %*%
      Bt_Yi %*% (t(invQ_PhiT_Yi) %*% diag(delta/lambda, q, q) %*% invQ)
    term3 <- 2 * sigma2 * BtB %*% Theta_invQ %*%
      diag(delta/lambda, q, q) %*% invQ

    H_Theta_eta <- H_Theta_eta +
      1 / n * (term1 + term2 + term3)

    ### H_eta_eta

    term1 <- delta / lambda * invQ_PhiT_Yi^2
    term2 <- -2 * sigma2 / lambda * invQ_PhiT_Yi *
      (invQ %*% (delta / lambda * invQ_PhiT_Yi))
    term3 <- -sigma2^2 / lambda * diag(invQ %*% diag(delta/lambda, q, q) %*% invQ)
    term4 <- sigma2 * diag(invQ) * delta / lambda

    H_eta_eta <- H_eta_eta + 1 / n * (term1 + term2 + term3 + term4)

    ### H_zeta_eta

    # term0 <- sum((invQ_PhiT_Yi * delta/lambda) * invQ_PhiT_Yi)
    term1 <- 2 * sigma2 * sum(
      (invQ_PhiT_Yi * delta/lambda) * invQ %*% (invQ_PhiT_Yi / lambda)
    )
    term2 <- sigma2^2 * sum(diag(invQ %*% diag(delta/lambda, q, q) %*% invQ) / lambda)
    term3 <- -sigma2 * sum(diag(invQ) * delta / lambda)

    H_zeta_eta <- H_zeta_eta + 1 / n * sum(term1 + term2 + term3)
  }

  list(
    H_Theta_eta = as.matrix(invG %*% H_Theta_eta),
    H_eta_eta = as.vector(H_eta_eta),
    H_zeta_eta = as.numeric(H_zeta_eta)
  )
}

delta_eta <- rnorm(q)
Hess_x_eta_apply(delta_eta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2)

check_H_Theta_eta <- check.derivatives(
  0,
  \(h) {
    compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda * exp(h * delta_eta), sigma2
    )[["grad_Theta"]] |> as.vector()
  },
  \(h) {
    Hess_x_eta_apply(
      delta_eta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_Theta_eta"]] |> as.vector()
  }
)

check_H_eta_eta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda * exp(h * delta_eta), sigma2
    )[["grad_eta"]] |> as.vector()
  },
  \(h) {
    Hess_x_eta_apply(
      delta_eta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_eta_eta"]] |> as.vector()
  }
)

check_H_zeta_eta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda * exp(h * delta_eta), sigma2
    )[["grad_zeta"]] |> as.vector()
  },
  \(h) {
    Hess_x_eta_apply(
      delta_eta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_zeta_eta"]] |> as.vector()
  }
)


# Check Hessian zeta -----------------------------------------------------

Hess_x_zeta_apply <- function(
  delta,
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(length(delta) == 1)
  stopifnot(length(Ltid) == n)
  H_Theta_zeta <- matrix(0, p, q)
  H_eta_zeta <- numeric(q)
  H_zeta_zeta <- 0

  for (i in seq_len(n)) {
    mi <- length(Ltid[[i]])
    Bi <- B[Ltid[[i]], , drop = F]
    Yi <- Ly[[i]]
    if (!is.null(theta_mu)) {
      Yi <- Yi - Bi %*% theta_mu
    }
    BtB <- as.matrix(crossprod(Bi))
    Phi <- as.matrix(Bi %*% Theta)
    Bt_Yi <- as.matrix(crossprod(Bi, Yi))
    PhiT_Yi <- crossprod(Phi, Yi)
    Phit_Phi <- crossprod(Phi)
    Q <- crossprod(Phi, Phi) + sigma2 * diag(1 / lambda, q, q)
    invQ <- solve(Q)
    Theta_invQ <- Theta %*% invQ
    invQ_PhiT_Yi <- invQ %*% PhiT_Yi

    ### H_Theta_zeta

    term1 <- -2 * delta / sigma2 * (
      (BtB %*% Theta_invQ %*% t(Theta) - diag(1, p, p)) %*%
        Bt_Yi %*% t(invQ_PhiT_Yi)
    )
    term2 <- -2 * delta * BtB %*% Theta_invQ %*% diag(1/lambda, q, q) %*%
      invQ_PhiT_Yi %*% t(invQ_PhiT_Yi)
    term3 <- -2 * delta * (BtB %*% Theta_invQ %*% t(Theta) - diag(1, p, p)) %*%
      Bt_Yi %*% (t(invQ_PhiT_Yi) %*% diag(1/lambda, q, q) %*% invQ)
    term4 <- -2 * sigma2 * delta * BtB %*% Theta_invQ %*%
      diag(1/lambda, q, q) %*% invQ

    H_Theta_zeta <- H_Theta_zeta +
      1 / n * (term1 + term2 + term3 + term4)

    ### H_eta_zeta

    term1 <- 2 * delta * sigma2 * (invQ_PhiT_Yi / lambda) * (invQ %*% (invQ_PhiT_Yi / lambda))
    term2 <- -sigma2 * delta * diag(invQ) / lambda
    term3 <- sigma2^2 * delta * diag(invQ %*% diag(1/lambda,q,q) %*% invQ) / lambda

    H_eta_zeta <- H_eta_zeta + 1 / n * (term1 + term2 + term3)

    ### H_zeta_zeta

    term1 <- delta / sigma2 * sum(Yi^2)
    term2 <- -delta / sigma2 * sum(PhiT_Yi * invQ_PhiT_Yi)
    term3 <- -delta * sum(invQ_PhiT_Yi^2 / lambda)
    term4 <- -2 * delta * sigma2 *
      sum(invQ_PhiT_Yi / lambda * invQ %*% (invQ_PhiT_Yi / lambda))
    term5 <- sigma2 * delta * sum(diag(invQ) / lambda)
    term6 <- -sigma2^2 * delta * sum(diag(invQ %*% diag(1/lambda,q,q) %*% invQ) / lambda)

    H_zeta_zeta <- H_zeta_zeta + 1 / n *
      sum(term1 + term2 + term3 + term4 + term5 + term6)
  }

  list(
    H_Theta_zeta = as.matrix(invG %*% H_Theta_zeta),
    H_eta_zeta = as.vector(H_eta_zeta),
    H_zeta_zeta = as.numeric(H_zeta_zeta)
  )
}

delta_zeta <- rnorm(1)
Hess_x_zeta_apply(delta_zeta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2)

check_H_Theta_zeta <- check.derivatives(
  0,
  \(h) {
    compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2 * exp(h * delta_zeta)
    )[["grad_Theta"]] |> as.vector()
  },
  \(h) {
    Hess_x_zeta_apply(
      delta_zeta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_Theta_zeta"]] |> as.vector()
  }
)

check_H_Theta_zeta <- check.derivatives(
  0,
  \(h) {
    as.vector(
      G %*% compute_grads(
        dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2 * exp(h * delta_zeta)
      )[["grad_Theta"]]
    )
  },
  \(h) {
    as.vector(
      G %*% Hess_x_zeta_apply(
        delta_zeta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
      )[["H_Theta_zeta"]]
    )
  }
)

check_H_eta_zeta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2 * exp(h * delta_zeta)
    )[["grad_eta"]] |> as.vector()
  },
  \(h) {
    Hess_x_zeta_apply(
      delta_zeta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_eta_zeta"]] |> as.vector()
  }
)

check_H_zeta_zeta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2 * exp(h * delta_zeta)
    )[["grad_zeta"]] |> as.vector()
  },
  \(h) {
    Hess_x_zeta_apply(
      delta_zeta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    )[["H_zeta_zeta"]] |> as.vector()
  }
)


# Vectorized Hessian -----------------------------------------------

get_commute_index <- function(p, q) {
  i <- rep(seq_len(p), times = q)   # row indices in X
  j <- rep(seq_len(q), each  = p)   # col indices in X
  ridx <- j + (i - 1L) * q             # positions in vec(t(X))
  cidx <- i + (j - 1L) * p             # positions in vec(X)
  perm <- integer(p * q)
  perm[ridx] <- cidx
  invperm <- integer(p * q)
  invperm[perm] <- seq_len(p * q)
  invperm
}

perm <- get_commute_index(p, q)


row_outer_prod <- function(A, B) {
  stopifnot(is.matrix(A) && is.matrix(B))
  stopifnot(nrow(A) == nrow(B))
  sapply(
    seq_len(nrow(A)), \(i) {
      c(outer(A[i,], B[i,]))
    }
  ) |> matrix(ncol = nrow(A)) |> t()
}

computeVecHess <- function(
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(ncol(Delta) == q)
  stopifnot(nrow(Delta) == p)
  stopifnot(length(Ltid) == n)
  H_Theta_Theta <- array(0, dim=c(p*q, p*q))
  H_eta_Theta <- array(0, dim=c(q, p*q))
  H_zeta_Theta <- array(0, dim=c(p*q))
  H_Theta_eta <- matrix(0, p*q, q)
  H_eta_eta <- matrix(0, q, q)
  H_zeta_eta <- numeric(q)
  H_Theta_zeta <- numeric(p*q)
  H_eta_zeta <- numeric(q)
  H_zeta_zeta <- numeric(1)

  for (i in seq_len(n)) {
    mi <- length(Ltid[[i]])
    Bi <- B[Ltid[[i]], , drop = F]
    Yi <- Ly[[i]]
    if (!is.null(theta_mu)) {
      Yi <- Yi - Bi %*% theta_mu
    }
    Phi <- as.matrix(Bi %*% Theta)
    Bt_Yi <- as.matrix(crossprod(Bi, Yi))
    BtB <- as.matrix(crossprod(Bi, Bi))
    PhiT_Yi <- crossprod(Phi, Yi)
    Phit_Phi <- crossprod(Phi)
    BtPhi <- as.matrix(crossprod(Bi, Phi))
    Q <- crossprod(Phi, Phi) + sigma2 * diag(1 / lambda, q, q)
    invQ <- solve(Q)
    invQinvLam <- invQ %*% diag(1 / lambda, q, q)
    Phi_invQ <- Phi %*% invQ
    invQ_PhiT_Yi <- invQ %*% PhiT_Yi
    tmp <- BtPhi %*% invQ_PhiT_Yi - Bt_Yi

    ### H_Theta_Theta
    fA <- 2 * tau * as.matrix(Omega)
    fC1 <- 2 / sigma2 * (BtB - BtPhi %*% invQ %*% t(BtPhi))
    fD1 <- invQ_PhiT_Yi %*% t(invQ_PhiT_Yi) + sigma2 * invQ
    fC2 <- -2 / sigma2 * tmp %*% t(tmp)
    fD2 <- invQ
    fE1 <- 2 / sigma2 * BtPhi %*% invQ
    fF1 <- -tmp %*% t(invQ_PhiT_Yi)
    fE2 <- -2 * (tmp %*% t(invQ_PhiT_Yi) / sigma2 + BtPhi %*% invQ)
    fF2 <- BtPhi %*% invQ

    H_Theta_Theta <- H_Theta_Theta + 1 / n * (
      diag(1,q,q) %x% as.matrix(invG %*% fA)
        + t(fD1) %x% as.matrix(invG %*% fC1)
        + t(fD2) %x% as.matrix(invG %*% fC2)
    )
    H_Theta_Theta <- H_Theta_Theta + 1 / n * (
        + t(fF1) %x% as.matrix(invG %*% fE1)
        + t(fF2) %x% as.matrix(invG %*% fE2)
    )[,perm]

    ### H_eta_Theta
    fC <- 2 * diag(c(invQ_PhiT_Yi) / lambda, q, q) %*% invQ %*% t(BtPhi)
    fD <- invQ_PhiT_Yi
    fE <- 2 * diag(c(invQ_PhiT_Yi) / lambda, q, q) %*% invQ
    fF <- BtPhi %*% invQ_PhiT_Yi - Bt_Yi
    fa <- invQ %*% t(BtPhi)
    fb <- invQ

    H_eta_Theta <- H_eta_Theta + 1 / n * (
      t(fD) %x% fC
      + 2 * sigma2 * diag(1/lambda,q,q) %*% row_outer_prod(fa, fb)
    )
    H_eta_Theta <- H_eta_Theta + 1 / n * (
      t(fF) %x% fE
    )[,perm]

    ### H_zeta_Theta
    fA <- -2 * sigma2 * BtPhi %*% invQ %*%
      diag(1/lambda,q,q) %*% invQ
    fC1 <- 2 * as.vector(Bt_Yi)
    fD1 <- as.vector(
      invQ_PhiT_Yi / sigma2 + invQ %*% (invQ_PhiT_Yi / lambda)
    )
    fC2 <- -as.vector(BtPhi %*% invQ_PhiT_Yi)
    fD2 <- as.vector(
      invQ_PhiT_Yi / sigma2
        + 2 * invQ %*% (invQ_PhiT_Yi / lambda)
    )
    fE <- -as.vector(invQ_PhiT_Yi)
    fF <- as.vector(
      BtPhi %*% (
      invQ_PhiT_Yi / sigma2
        + 2 * invQ %*% (invQ_PhiT_Yi / lambda)
    ))

    H_zeta_Theta <- H_zeta_Theta + 1 / n * (
      as.vector(fA)
      + fD1 %x% fC1
      + fD2 %x% fC2
    )
    H_zeta_Theta <- H_zeta_Theta + 1 / n * (
      + fF %x% fE
    )[perm]

    ### H_Theta_eta
    fA1 <- 2 * invG %*% BtPhi %*% invQinvLam
    fB1 <- (invQ_PhiT_Yi %*% t(invQ_PhiT_Yi) + sigma2 * invQ)
    fA2 <- 2 * invG %*% (BtPhi %*% invQ %*% t(Theta) - diag(1,p,p)) %*%
      Bt_Yi %*% t(invQ_PhiT_Yi/lambda)
    fB2 <- invQ

    for (k in seq_len(q)) {
      tmpA <- cbind(fA1[,k], fA2[,k])
      tmpB <- rbind(fB1[k,], fB2[k,])
      H_Theta_eta[,k] <- H_Theta_eta[,k] + 1 / n *
        as.vector(tmpA %*% tmpB)
    }

    ### H_eta_eta
    fA1 <- diag(c(invQ_PhiT_Yi) / lambda, q, q) %*%
      (diag(1,q,q) - 2 * sigma2 * invQ %*% diag(1/lambda,q,q)) %*%
      diag(c(invQ_PhiT_Yi), q, q)
    fA2 <- -sigma2^2 * invQinvLam * t(invQinvLam) +
      sigma2 * invQinvLam * diag(1,q,q)
    
    H_eta_eta <- H_eta_eta + 1 / n * (fA1 + fA2)

    ### H_zeta_eta
    fa1 <- 2 * sigma2 / lambda * as.vector(
      invQ_PhiT_Yi * invQinvLam %*% invQ_PhiT_Yi
    )
    fa2 <- sigma2^2 * as.vector(diag(invQinvLam %*% invQinvLam))
    fa3 <- -sigma2 * as.vector(diag(invQinvLam))

    H_zeta_eta <- H_zeta_eta + 1 / n * (fa1 + fa2 + fa3)

    ### H_Theta_zeta
    H_Theta_zeta <- H_Theta_zeta + 1 / n * as.vector(invG %*% (
      -2 / sigma2 * tmp %*% t(invQ_PhiT_Yi)
      -2 * BtPhi %*% invQinvLam %*% invQ_PhiT_Yi %*% t(invQ_PhiT_Yi)
      -2 * tmp %*% t(invQinvLam %*% invQ_PhiT_Yi)
      -2 * sigma2 * BtPhi %*% invQinvLam %*% invQ
    ))

    ### H_eta_zeta
    H_eta_zeta <- H_eta_zeta + 1 / n * as.vector(
      2 * sigma2 / lambda *
        invQ_PhiT_Yi * invQinvLam %*% invQ_PhiT_Yi
        -sigma2 * diag(invQinvLam)
        +sigma2^2 * diag(invQinvLam %*% invQinvLam)
    )

    ### H_zeta_zeta
    H_zeta_zeta <- H_zeta_zeta + 1 / n * as.numeric(
      sum(Yi^2) / sigma2
        - sum(PhiT_Yi * invQ_PhiT_Yi) / sigma2
        - sum(invQ_PhiT_Yi^2 / lambda)
        - 2 * sigma2 * sum((invQ_PhiT_Yi / lambda) * invQinvLam %*% invQ_PhiT_Yi)
        + sigma2 * sum(diag(invQinvLam))
        - sigma2^2 * sum(diag(invQinvLam %*% invQinvLam))
    )
  }

  out <- list(
    H_Theta_Theta = H_Theta_Theta,
    H_eta_Theta = H_eta_Theta,
    H_zeta_Theta = H_zeta_Theta,
    H_Theta_eta = H_Theta_eta,
    H_eta_eta = H_eta_eta,
    H_zeta_eta = H_zeta_eta,
    H_Theta_zeta = H_Theta_zeta,
    H_eta_zeta = H_eta_zeta,
    H_zeta_zeta = H_zeta_zeta
  )
  class(out) <- "vecHess"
  out
}

vecHess_apply <- function(
  vh, Delta, delta_eta, delta_zeta,
  grad_Theta = NULL,
  manifold = TRUE
) {
  stopifnot(is(vh, "vecHess"))
  if (manifold) {
    stopifnot(!is.null(grad_Theta))
    # browser()
    H_Th_Th_Del <- matrix(
      vh[['H_Theta_Theta']] %*% as.vector(Delta),
      p, q
    )
    symThGgTh <- sym(Theta, grad_Theta, G)
    symDelGgTh <- sym(Delta, grad_Theta, G)
    H_Th_Th_Del <- manifold.Stiefel.project(H_Th_Th_Del, Theta, G) -
      Delta %*% symThGgTh - Theta %*% symDelGgTh
    H_Th_Th_Del <- manifold.Stiefel.project(H_Th_Th_Del, Theta, G)

    H_Th_eta_del <- matrix(vh[['H_Theta_eta']] %*% delta_eta, p, q)
    H_Th_eta_del <- manifold.Stiefel.project(H_Th_eta_del, Theta, G)

    H_Th_zeta_del <- matrix(vh[['H_Theta_zeta']] * delta_zeta, p, q)
    H_Th_zeta_del <- manifold.Stiefel.project(H_Th_zeta_del, Theta, G)

    H_Th_Psi_Del <- H_Th_Th_Del + H_Th_eta_del + H_Th_zeta_del
  } else {
    H_Th_Psi_Del <- matrix(
      vh[['H_Theta_Theta']] %*% as.vector(Delta)
        + vh[['H_Theta_eta']] %*% delta_eta
        + vh[['H_Theta_zeta']] * delta_zeta,
      p, q
    )
  }
  H_eta_Psi_Del <- as.vector(
    vh[['H_eta_Theta']] %*% as.vector(Delta)
    + vh[['H_eta_eta']] %*% delta_eta
    + vh[['H_eta_zeta']] * delta_zeta
  )
  H_zeta_Psi_Del <- (
    sum(vh[['H_zeta_Theta']] * as.vector(Delta))
    + sum(vh[['H_zeta_eta']] * delta_eta)
    + vh[['H_zeta_zeta']] * delta_zeta
  )
  out <- list(
    H_Theta_Psi_Delta = H_Th_Psi_Del,
    H_eta_Psi_Delta = H_eta_Psi_Del,
    H_zeta_Psi_Delta = H_zeta_Psi_Del
  )
  if (manifold) {
    class(out) <- c("vecHessApply", "manifold")
  } else {
    class(out) <- c("vecHessApply", "ambient")
  }
  out
}


Delta <- matrix(rnorm(p*q), p, q)
Delta <- manifold.Stiefel.project(Delta, Theta, G)
delta_eta <- rnorm(q)
delta_zeta <- rnorm(1)
grad_Theta <- compute_grads(
  dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
)[["grad_Theta"]]
grad_eta <- compute_grads(
  dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
)[["grad_eta"]]
grad_zeta <- objfun(
  dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
)[["grad_zeta"]]


vh <- computeVecHess(
  dat$Ly[1:N0], dat$Ltid[1:N0],
  Theta, lambda, sigma2
)
vha <- vecHess_apply(vh, Delta, delta_eta, delta_zeta, grad_Theta, FALSE)


H_x_th_Del <- Hess_x_Theta_apply(
  Delta,
  dat$Ly[1:N0], dat$Ltid[1:N0],
  Theta, lambda, sigma2
)
H_x_eta_del <- Hess_x_eta_apply(
  delta_eta,
  dat$Ly[1:N0], dat$Ltid[1:N0],
  Theta, lambda, sigma2
)
H_x_zeta_del <- Hess_x_zeta_apply(
  delta_zeta,
  dat$Ly[1:N0], dat$Ltid[1:N0],
  Theta, lambda, sigma2
)
H_x_th_Del$H_Theta_Theta + H_x_eta_del$H_Theta_eta + H_x_zeta_del$H_Theta_zeta
vha$H_Theta_Psi_Delta
H_x_th_Del$H_eta_Theta + H_x_eta_del$H_eta_eta + H_x_zeta_del$H_eta_zeta
vha$H_eta_Psi_Delta
H_x_th_Del$H_zeta_Theta + H_x_eta_del$H_zeta_eta + H_x_zeta_del$H_zeta_zeta
vha$H_zeta_Psi_Delta


# Riemannian Hessian -----------------------------------------------------

Hess_Theta_apply <- function(
  Delta,
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(ncol(Delta) == q)
  stopifnot(nrow(Delta) == p)
  stopifnot(length(Ltid) == n)
  Hess_Theta <- matrix(0, p, q)

  for (i in seq_len(n)) {
    mi <- length(Ltid[[i]])
    Bi <- B[Ltid[[i]], , drop = F]
    Yi <- Ly[[i]]
    if (!is.null(theta_mu)) {
      Yi <- Yi - Bi %*% theta_mu
    }
    Phi <- as.matrix(Bi %*% Theta)
    Bt_Yi <- as.matrix(crossprod(Bi, Yi))
    BtPhi <- as.matrix(crossprod(Bi, Phi))
    PhiT_Yi <- crossprod(Phi, Yi)
    Phit_Phi <- crossprod(Phi)
    Q <- crossprod(Phi, Phi) + sigma2 * diag(1 / lambda, q, q)
    invQ <- solve(Q)
    Theta_invQ <- Theta %*% invQ
    invQ_PhiT_Yi <- invQ %*% PhiT_Yi
    BDel <- as.matrix(Bi %*% Delta)
    BDelT_Yi <- crossprod(BDel, Yi)
    Phit_BDel <- crossprod(Phi, BDel)
    BtB <- as.matrix(crossprod(Bi))
    Delta_invQ <- Delta %*% invQ
    invQ_BDelT_Yi <- invQ %*% BDelT_Yi

    factor1 <- BtPhi %*% t(Theta_invQ) - diag(1,p,p)
    factor2 <- Bt_Yi %*% t(invQ_PhiT_Yi)

    term1 <- 2 / sigma2 * BtB %*% (
      Delta_invQ %*% t(Theta) + Theta_invQ %*% t(Delta)
      - Theta_invQ %*% (Phit_BDel + t(Phit_BDel)) %*% t(Theta_invQ) 
    ) %*% factor2
    term2 <- 2 / sigma2 * factor1 %*%
      Bt_Yi %*% (t(Bt_Yi) %*% (
        Delta - Theta_invQ %*% (Phit_BDel + t(Phit_BDel))
      ) %*% invQ)
    term3 <- 2 * BtB %*% {
      Delta - Theta_invQ %*% (Phit_BDel + t(Phit_BDel))
    } %*% invQ
    term4 <- 2 * tau * as.matrix(Omega) %*% Delta

    Hess_Theta <- Hess_Theta +
      1 / n * (term1 + term2 + term3 + term4)
  }

  as.matrix(invG %*% Hess_Theta)
}


maniHess_Theta_apply <- function(
  Delta,
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3
) {
  Hess_Theta <- Hess_Theta_apply(
    Delta, Ly, Ltid,
    Theta, lambda, sigma2, theta_mu, tau
  )
  grad_Theta <- compute_grads(
    Ly, Ltid, Theta, lambda, sigma2, theta_mu, tau
  )[['grad_Theta']]
  ThGZ <- as.matrix(t(Theta) %*% G %*% grad_Theta)
  DelGZ <- as.matrix(t(Delta) %*% G %*% grad_Theta)
  Dproj <- - 0.5 * Delta %*% (ThGZ + t(ThGZ)) -
    0.5 * Theta %*% (DelGZ + t(DelGZ))

  manifold.Stiefel.project(Hess_Theta, Theta, G) + Dproj
}

Delta <- matrix(rnorm(p*q), p, q)
Delta <- manifold.Stiefel.project(Delta, Theta, G)
Hess_Theta_apply(Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2)
maniHess_Theta_apply(Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2)

check_H_Theta <- check.derivatives(
  0,
  \(h) {
    grad_Theta <- compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta + h * Delta, lambda, sigma2
    )[["grad_Theta"]]
    grad_Theta
  },
  \(h) {
    as.vector(Hess_Theta_apply(
      Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    ))
  }
)

grad_Theta <- compute_grads(
  dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
)[["grad_Theta"]]

check_grad_proj <- check.derivatives(
  0,
  \(h) {
    manifold.Stiefel.project(grad_Theta, Theta + h * Delta, G) |> as.vector()
  },
  \(h) {
    Z <- grad_Theta
    ThGZ <- as.matrix(t(Theta) %*% G %*% grad_Theta)
    DelGZ <- as.matrix(t(Delta) %*% G %*% grad_Theta)
    Dproj <- - 0.5 * Delta %*% (ThGZ + t(ThGZ)) -
      0.5 * Theta %*% (DelGZ + t(DelGZ))
    Dproj |> as.vector()
  }
)

check_maniH_Theta <- check.derivatives(
  0,
  \(h) {
    grad_Theta <- compute_grads(
      dat$Ly[1:N0], dat$Ltid[1:N0], Theta + h * Delta, lambda, sigma2
    )[["grad_Theta"]]
    as.vector(manifold.Stiefel.project(grad_Theta, Theta + h * Delta, G))
  },
  \(h) {
    maniHess_Theta_apply(
      Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
    ) |> as.vector()
  }
)

vh <- computeVecHess(
  dat$Ly[1:N0], dat$Ltid[1:N0],
  Theta, lambda, sigma2
)
vha <- vecHess_apply(vh, Delta, delta_eta, delta_zeta, grad_Theta, TRUE)
vha$H_Theta_Psi_Delta
RHessTh <- maniHess_Theta_apply(
  Delta, dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
) |> manifold.Stiefel.project(Theta, G)
RHessTh +
  manifold.Stiefel.project(H_x_eta_del$H_Theta_eta, Theta, G) +
  manifold.Stiefel.project(H_x_zeta_del$H_Theta_zeta, Theta, G)


# Cov of grad ------------------------------------------------------------
cov_grad_Theta <- foreach(
  i = iterators::iter(1:1000),
  .combine = "+"
) %do% {
  rgrad_Theta <- compute_grads(
    dat$Ly[i], dat$Ltid[i], Theta, lambda, sigma2
  )[["grad_Theta"]] |> 
    manifold.Stiefel.project(Theta, G)
  outer(c(rgrad_Theta), c(rgrad_Theta))
} / 1000
image(Matrix::Matrix(cov_grad_Theta))


# Check PSD of Riemannian Hessian ----------------------------------------

grad_Theta <- compute_grads(
  dat$Ly[1:N0], dat$Ltid[1:N0], Theta, lambda, sigma2
)[["grad_Theta"]]
vh <- computeVecHess(
  dat$Ly[1:N0], dat$Ltid[1:N0],
  Theta, lambda, sigma2
)

Delta <- matrix(rnorm(p*q), p, q)
Delta <- manifold.Stiefel.project(Delta, Theta, G)
delta_eta <- rnorm(q)
delta_zeta <- rnorm(1)

Delta2 <- matrix(rnorm(p*q), p, q)
Delta2 <- manifold.Stiefel.project(Delta2, Theta, G)
delta2_eta <- rnorm(q)
delta2_zeta <- rnorm(1)

vha <- vecHess_apply(vh, Delta, delta_eta, delta_zeta, grad_Theta, TRUE)
vha2 <- vecHess_apply(vh, Delta2, delta2_eta, delta2_zeta, grad_Theta, TRUE)

sum(G %*% vha$H_Theta_Psi_Delta * Delta) +
  sum(vha$H_eta_Psi_Delta * delta_eta) +
  sum(vha$H_zeta_Psi_Delta * delta_zeta)
sum(G %*% vha$H_Theta_Psi_Delta * Delta2) +
  sum(vha$H_eta_Psi_Delta * delta2_eta) +
  sum(vha$H_zeta_Psi_Delta * delta2_zeta)
sum(G %*% vha2$H_Theta_Psi_Delta * Delta) +
  sum(vha2$H_eta_Psi_Delta * delta_eta) +
  sum(vha2$H_zeta_Psi_Delta * delta_zeta)


# Manually implemented CG ------------------------------------------------

# Conjugate Gradient for linear operator equations H(x) = b
# Assumptions: H is self-adjoint (symmetric) and positive definite
# Interface:
#   H        : function(x) returning H(x)
#   b        : RHS vector (same shape as x)
#   x0       : initial guess (defaults to zeros like b)
#   tol      : relative residual tolerance ||r||/||b||
#   maxit    : max iterations
#   Msolve   : preconditioner solve function z = M^{-1} r  (identity if NULL)
#   ip       : inner product function <u,v> (defaults to dot product)
#   callback : optional function(it, x, rnorm_rel) called each iter
# Returns: list(x=..., converged=..., iter=..., relres=..., resnorm_history=...)

cg_operator <- function(H, b, x0 = NULL, tol = 1e-8, maxit = NULL,
                        Msolve = NULL,
                        ip = function(u, v) drop(sum(u * v)),
                        callback = NULL) {
  # shape helpers
  zlike <- function(v) { v[seq_along(b)] }  # keep shape like b
  nrm2  <- function(v) sqrt(max(ip(v, v), 0))

  if (is.null(x0)) x <- zlike(0*b) else x <- zlike(x0)
  if (is.null(maxit)) maxit <- 2L * length(b)  # generous default for safety
  if (is.null(Msolve)) Msolve <- function(r) r  # identity preconditioner

  r <- b - H(x)                 # initial residual
  bnorm <- nrm2(b)
  if (bnorm == 0) bnorm <- 1.0  # avoid div by zero if b=0
  relres <- nrm2(r) / bnorm

  # Early exit
  if (relres <= tol) {
    return(list(x = x, converged = TRUE, iter = 0L,
                relres = relres, resnorm_history = relres))
  }

  z  <- Msolve(r)
  p  <- z
  rz_old <- ip(r, z)

  if (!is.finite(rz_old) || rz_old <= 0) {
    stop("Preconditioner or inner product produced non-positive <r, M^{-1}r> at start.")
  }

  history <- numeric()
  history[1] <- relres
  converged <- FALSE
  it <- 0L

  for (k in 1:maxit) {
    Hp <- H(p)
    pHp <- ip(p, Hp)

    # SPD safeguard
    if (!is.finite(pHp) || pHp <= 0) {
      warning(sprintf("CG breakdown (non-positive curvature) at iter %d: <p, H p> = %g", k, pHp))
      break
    }

    alpha <- rz_old / pHp
    x <- x + alpha * p
    r <- r - alpha * Hp

    relres <- nrm2(r) / bnorm
    history[k + 1] <- relres
    it <- k
    if (!is.null(callback)) try(callback(it, x, relres), silent = TRUE)

    if (relres <= tol) { converged <- TRUE; break }

    z <- Msolve(r)
    rz_new <- ip(r, z)

    if (!is.finite(rz_new) || rz_new <= 0) {
      warning(sprintf("CG breakdown at iter %d: <r, M^{-1} r> = %g", k, rz_new))
      break
    }

    beta <- rz_new / rz_old
    p <- z + beta * p
    rz_old <- rz_new
  }

  list(x = x, converged = converged, iter = it,
       relres = relres, resnorm_history = history)
}


vh <- computeVecHess(
  dat$Ly[1:N], dat$Ltid[1:N],
  Theta, lambda, sigma2
)
grad_Theta <- compute_grads(
  dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
)[["grad_Theta"]]
i <- sample(1:N0, 1)
rgrad_i_Theta <- compute_grads(
  dat$Ly[i], dat$Ltid[i], Theta, lambda, sigma2
)[["grad_Theta"]] |> 
  manifold.Stiefel.project(Theta, G)

vecHess_apply_wrapper <- function(x) {
  x <- vecProj(x)
  Del_Th <- matrix(x[1:(p*q)], p ,q)
  del_eta <- x[(p*q+1):(p*q+q)]
  del_zeta <- x[(p*q+q+1):(p*q+q+1)]
  vha <- vecHess_apply(vh, Del_Th, del_eta, del_zeta, grad_Theta, TRUE)
  return(c(
    vha$H_Theta_Psi_Delta,
    vha$H_eta_Psi_Delta,
    vha$H_zeta_Psi_Delta
  ))
}

vecProj <- function(x) {
  Delta <- matrix(x[1:(p*q)], p, q)
  Delta <- manifold.Stiefel.project(Delta, Theta, G)
  c(Delta, tail(x, q + 1))
}

vecGrad <- c(rgrad_i_Theta, grad_eta, grad_zeta)

# FIXME: CG algo still does not work. Fix it.

cgres <- cg_operator(
  vecHess_apply_wrapper,
  b = vecGrad,
  x0 = vecGrad
)

vecHess_apply_wrapper(vecProj(cgres$x))
vecGrad

# An alternative: solve the ls problem with an nls algo.

opres <- gslnls::gsl_nls(
  fn = \(x) { vecHess_apply_wrapper(vecProj(x)) },
  y = vecGrad,
  start = vecGrad
)

vecU <- vecProj(opres$m$getPars())
vecHess_apply_wrapper(vecU)
vecGrad


# Differential of QR -----------------------------------------------------

diff_Qfactor <- function(dA, A) {
  # differential of the Q factor for square matrix A
  p <- nrow(A)
  qrObj <- qr(A)
  Q <- qr.Q(qrObj)
  R <- qr.R(qrObj)
  invR <- solve(R)
  QtdAinvR <- crossprod(Q, dA %*% invR)
  L <- QtdAinvR
  L[upper.tri(L, diag = TRUE)] <- 0
  (diag(1,p,p) - Q %*% t(Q)) %*% dA %*% invR + Q %*% (L - t(L))
}

A <- matrix(rnorm(p*q), p, q)
dA <- matrix(rnorm(p*q), p, q)
diff_Qfactor(A, dA)

check_diff_Q <- check.derivatives(
  0, \(h) {qr.Q(qr(A + h * dA)) |> as.vector()},
  \(h) {diff_Qfactor(dA, A) |> as.vector()}
)

V2phi <- function(V, Theta) {
  # Theta <- manifold.Stiefel.retract(Theta + V, G)
  Theta <- sqrtGinv %*% qr.Q(qr(sqrtG %*% (Theta + V)))
  as.matrix(B %*% Theta)
}

dV2phi <- function(dV, V, Theta) {
  A <- sqrtG %*% (Theta + V)
  dA <- sqrtG %*% dV
  dQ <- diff_Qfactor(dA, A)
  as.matrix(B %*% sqrtGinv %*% dQ)
}

dV2phi.test <- function(dV, V, Theta) {
  p <- nrow(Theta)
  A <- sqrtG %*% (Theta + V)
  dA <- sqrtG %*% dV
  qrObj <- qr(A)
  Q <- qr.Q(qrObj)
  R <- qr.R(qrObj)
  invR <- solve(R)
  QtdAinvR <- crossprod(Q, dA %*% invR)
  L <- QtdAinvR
  L[upper.tri(L, diag = TRUE)] <- 0
  dQ <- (
    (diag(1,p,p) - Q %*% t(Q)) %*% dA %*% invR
      + Q %*% (L - t(L))
  )
  as.matrix(B %*% sqrtGinv %*% dQ)
}

d0V2phi.test <- function(dV, Theta) {
  p <- nrow(Theta)
  q <- ncol(Theta)
  Q <- sqrtG %*% Theta
  R <- diag(1, q, q)
  invR <- R
  dA <- sqrtG %*% dV
  QtdAinvR <- crossprod(Q, dA %*% invR)
  L <- QtdAinvR
  L[upper.tri(L, diag = TRUE)] <- 0
  dQ <- (
    (diag(1,p,p) - Q %*% t(Q)) %*% dA %*% invR
      + Q %*% (L - t(L))
  )
  as.matrix(B %*% sqrtGinv %*% dQ)
}

grad0_V2phi <- function(Theta) {
  Q <- sqrtG %*% Theta
  R <- diag(1, q, q)
  alphas <- as.matrix(B %*% sqrtGinv) # m x p, rows for alpha
  gammas <- as.matrix(B %*% Theta) # m x q, rows for gamma, =Phi
  ltri.idx <- lower.tri(matrix(0,q,q), diag=FALSE)
  utri.idx <- upper.tri(matrix(0,q,q), diag=FALSE)
  # S_{kt}: (eval, gamma, ek, k)
  Skt_array <- array(0, dim = c(nrow(B), q, q, q))
  # g_V: (eval, basis, ek, k)
  gVarr <- array(0, dim = c(nrow(B), p, q, q))
  for (k in seq_len(q)) {
    gam_ek <- array(0, dim = c(nrow(B), q, q))
    gam_ek[,,k] <- gammas
    Skt_array[,,,k] <- apply(
      gam_ek, 1, \(x) {
        gam_ek_j_l <- gam_ek_j_u <- matrix(0, q, q)
        gam_ek_j_l[ltri.idx] <- x[ltri.idx]
        gam_ek_j_u[utri.idx] <- x[utri.idx]
        gam_ek_j_l - t(gam_ek_j_u)
      }
    ) |> array(dim = c(q, q, nrow(B))) |> # (gamma, ek, eval)
      aperm(c(3,1,2))
    gVarr[,,k,k] = gVarr[,,k,k] +
      as.matrix(B) %*% invG - as.matrix(B %*% Theta) %*% t(Theta)
    gVarr[,,,k] = gVarr[,,,k] + 
      apply(Skt_array[,,,k], 3, \(X) X %*% t(Theta)) |> # (eval * basis, ek)
      array(dim = c(nrow(B), p, q))
  }
  gVarr
}

V <- matrix(rnorm(p*q), p, q)
V <- manifold.Stiefel.project(V, Theta, G)
dV <- matrix(rnorm(p*q), p, q)
dV <- manifold.Stiefel.project(dV, Theta, G)

check_diff_V2phi <- check.derivatives(
  0, \(h) {V2phi(V + h * dV, Theta) |> as.vector()},
  \(h) {dV2phi(dV, V, Theta) |> as.vector()}
)

check_diff_V2phi <- check.derivatives(
  0, \(h) {V2phi(matrix(0,p,q) + h * dV, Theta) |> as.vector()},
  \(h) {dV2phi(dV, matrix(0,p,q), Theta) |> as.vector()}
)

# DgV <- d0V2phi.test(dV, Theta)
DgV <- dV2phi(dV, matrix(0,p,q), Theta)
gV <- grad0_V2phi(Theta)
gV_G_dV <- apply(gV, c(1,4), \(X) sum(X * as.matrix(G %*% dV)))


# Solve J Hess^{-1} \Sigma Hess^{-1} J^* ------------------------------------------

Theta <- ThetaTrue[,1:q]
lambda <- lambdaTrue[1:q]
sigma2 <- sigma2true

vh <- computeVecHess(
  dat$Ly, dat$Ltid,
  Theta, lambda, sigma2
)
grad_Theta <- compute_grads(
  dat$Ly, dat$Ltid, Theta, lambda, sigma2
)[["grad_Theta"]]
grad_eta <- compute_grads(
  dat$Ly, dat$Ltid, Theta, lambda, sigma2
)[["grad_eta"]]
grad_zeta <- compute_grads(
  dat$Ly, dat$Ltid, Theta, lambda, sigma2
)[["grad_zeta"]]

vecProj <- function(x) {
  Delta <- matrix(x[1:(p*q)], p, q)
  Delta <- manifold.Stiefel.project(Delta, Theta, G)
  c(Delta, tail(x, q + 1))
}

library(doFuture)
plan(multisession(workers = 8))

foreach(
  i = 1:1000,
  .combine = "+"
) %dofuture% {
  gradObj <- compute_grads(
    dat$Ly[i], dat$Ltid[i], Theta, lambda, sigma2
  )
  rgrad_Theta <- gradObj[["grad_Theta"]] |> 
    manifold.Stiefel.project(Theta, G)
  grad_eta <- gradObj[["grad_eta"]]
  grad_zeta <- gradObj[["grad_zeta"]]
  vecGrad <- c(rgrad_Theta, grad_eta, grad_zeta)
  opres <- gslnls::gsl_nls(
    fn = \(x) { 
      x <- vecProj(x)
      Del_Th <- matrix(x[1:(p*q)], p ,q)
      del_eta <- x[(p*q+1):(p*q+q)]
      del_zeta <- x[(p*q+q+1):(p*q+q+1)]
      vha <- vecHess_apply(vh, Del_Th, del_eta, del_zeta, grad_Theta, TRUE)
      return(c(
        vha$H_Theta_Psi_Delta,
        vha$H_eta_Psi_Delta,
        vha$H_zeta_Psi_Delta
      ))
    },
    y = vecGrad,
    start = vecGrad
  )
  vecU <- vecProj(opres$m$getPars())
  U <- matrix(vecU[1:(p*q)], p, q)
  JU <- dV2phi(U, matrix(0,p,q), Theta)
  JU^2
}

phisig2 <- matrix(c(
  8592.774, 7287.720, 6110.106, 5062.715, 4146.464,
  3360.294, 2701.128, 2163.924, 1741.806, 1426.289,
  1207.578, 1074.960, 1017.278, 1023.532, 1084.234,
  1191.141, 1335.805, 1509.180, 1701.414, 1901.800,
  2098.873, 2280.660, 2435.090, 2550.546, 2616.571,
  2624.733, 2571.366, 2463.239, 2311.371, 2128.256,
  1927.027, 1720.790, 1522.119, 1342.717, 1193.248,
  1083.328, 1021.685, 1016.484, 1075.906, 1210.739,
  1434.521, 1760.909, 2202.846, 2771.921, 3477.828,
  4327.931, 5326.932, 6476.638, 7775.835, 9220.267,
  10802.708,

  2448.346, 2131.123, 1856.625, 1625.692, 1438.110,
  1292.578, 1186.727, 1117.187, 1079.699, 1069.284,
  1080.453, 1107.468, 1144.658, 1186.813, 1230.475,
  1274.075, 1316.767, 1358.052, 1397.527, 1434.683,
  1468.785, 1498.806, 1523.436, 1541.153, 1550.361,
  1549.597, 1538.307, 1518.477, 1492.752, 1463.585,
  1432.992, 1402.390, 1372.525, 1343.476, 1314.757,
  1285.492, 1254.684, 1221.565, 1186.084, 1150.589,
  1119.832, 1099.346, 1094.974, 1112.521, 1157.442,
  1234.596, 1348.046, 1500.918, 1695.303, 1932.221,
  2211.634,

  3123.213, 2626.612, 2281.931, 2057.748, 1927.797,
  1870.007, 1865.677, 1898.794, 1955.487, 2023.633,
  2092.591, 2153.088, 2197.240, 2218.815, 2215.774,
  2190.533, 2147.221, 2090.835, 2026.624, 1959.580,
  1894.057, 1833.516, 1780.381, 1736.024, 1700.872,
  1674.629, 1657.112, 1649.923, 1654.636, 1671.939,
  1701.393, 1741.341, 1788.929, 1840.271, 1890.739,
  1935.387, 1969.502, 1989.292, 1992.770, 1982.038,
  1963.024, 1943.417, 1932.365, 1940.398, 1979.491,
  2063.287, 2207.464, 2430.246, 2753.065, 3201.380,
  3805.630
), ncol = 3, byrow = FALSE) / 1000



par(mfrow=c(1,3))
for (k in seq_len(q)) {
  phikeval <- as.vector(B %*% Theta[,k])
  phik_u <- phikeval + 1.96 * sqrt(phisig2[,k]) 
  phik_l <- phikeval - 1.96 * sqrt(phisig2[,k]) 
  ylim = range(c( phik_u, phik_l))
  plot(tgrid, eval_fd(tgrid, phiTrueFunc)[,k], type="l",
    xlab="t", ylab=bquote(phi[.(k)]),
    ylim=ylim)
  lines(tgrid, phikeval, col=4)
  lines(tgrid, phik_u, col=2, lty=2)
  lines(tgrid, phik_l, col=2, lty=2)
}
par(mfrow=c(1,1))







# Solve Hessian operator -------------------------------------------------

