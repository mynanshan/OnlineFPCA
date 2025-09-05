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
invG <- sqrtGinv %*% sqrtGinv

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


# Check Derivatives ------------------------------------------------------

Theta = ThetaInit
lambda = lambdaInit
sigma2 = sigma2Init

library(nloptr)

N <- 20

check_grad_Theta <- check.derivatives(
  inits$Theta,
  \(x) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], matrix(x, p, q), inits$lambda, inits$sigma2,
      stats = "loss"
    )[["fval"]]
  },
  \(x) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], matrix(x, p, q), inits$lambda, inits$sigma2,
      stats = "grad"
    )[["grad_Theta"]] |> as.vector()
  }
)


check_grad_eta <- check.derivatives(
  log(inits$lambda),
  \(x) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], inits$Theta, exp(x), inits$sigma2,
      stats = "loss"
    )[["fval"]]
  },
  \(x) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], inits$Theta, exp(x), inits$sigma2,
      stats = "grad"
    )[["grad_eta"]]
  }
)


check_grad_zeta <- check.derivatives(
  log(inits$sigma2),
  \(x) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], inits$Theta, inits$lambda, exp(x),
      stats = "loss"
    )[["fval"]]
  },
  \(x) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], inits$Theta, inits$lambda, exp(x),
      stats = "grad"
    )[["grad_zeta"]]
  }
)


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

    term1 <- 2 / sigma2 * BtB %*% {
      Delta_invQ %*% t(Theta) + Theta_invQ %*% t(Delta)
      - Theta_invQ %*% (Phit_BDel + t(Phit_BDel)) %*% t(Theta_invQ) 
    } %*% Bt_Yi %*% (t(PhiT_Yi) %*% invQ)
    term2 <- 2 / sigma2 * (BtB %*% Theta_invQ %*% t(Theta) - diag(1,p,p)) %*%
      Bt_Yi %*% (t(Bt_Yi) %*% {
        Delta_invQ - Theta_invQ %*% (Phit_BDel + t(Phit_BDel)) %*% invQ
      })
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
    H_Theta_Theta = as.matrix(H_Theta_Theta),
    H_eta_Theta = as.vector(H_eta_Theta),
    H_zeta_Theta = as.numeric(H_zeta_Theta)
  )
}


Delta <- matrix(rnorm(p*q), p, q)
Delta <- manifold.Stiefel.project(Delta, Theta, G)
Hess_x_Theta_apply(Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2)

check_H_Theta_Theta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta + h * Delta, lambda, sigma2
    )[["grad_Theta"]] |> as.vector()
  },
  \(h) {
    Hess_x_Theta_apply(
      Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    )[["H_Theta_Theta"]] |> as.vector()
  }
)

check_H_eta_Theta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta + h * Delta, lambda, sigma2
    )[["grad_eta"]]
  },
  \(h) {
    Hess_x_Theta_apply(
      Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    )[["H_eta_Theta"]]
  }
)

check_H_zeta_Theta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta + h * Delta, lambda, sigma2
    )[["grad_zeta"]]
  },
  \(h) {
    Hess_x_Theta_apply(
      Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    )[["H_zeta_Theta"]]
  }
)


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

    H_Theta_eta <- H_Theta_eta + 1 / n * (term1 + term2 + term3)

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
    H_Theta_eta = H_Theta_eta,
    H_eta_eta = H_eta_eta,
    H_zeta_eta = H_zeta_eta
  )
}

delta_eta <- rnorm(q)
Hess_x_eta_apply(delta_eta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2)

check_H_Theta_eta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda * exp(h * delta_eta), sigma2
    )[["grad_Theta"]] |> as.vector()
  },
  \(h) {
    Hess_x_eta_apply(
      delta_eta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    )[["H_Theta_eta"]] |> as.vector()
  }
)

check_H_eta_eta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda * exp(h * delta_eta), sigma2
    )[["grad_eta"]] |> as.vector()
  },
  \(h) {
    Hess_x_eta_apply(
      delta_eta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    )[["H_eta_eta"]] |> as.vector()
  }
)

check_H_zeta_eta <- check.derivatives(
  0,
  \(h) {
    objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda * exp(h * delta_eta), sigma2
    )[["grad_zeta"]] |> as.vector()
  },
  \(h) {
    Hess_x_eta_apply(
      delta_eta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    )[["H_zeta_eta"]] |> as.vector()
  }
)


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


# Solve Hessian operator -------------------------------------------------


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

    factor1 <- crossprod(Bi, Phi) %*% t(Theta_invQ) - diag(1,p,p)
    factor2 <- crossprod(Bi, Yi) %*% t(invQ_PhiT_Yi)

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
    # term1 <- 0
    # term2 <- 0
    # term3 <- 0
    # term4 <- 0

    Hess_Theta <- Hess_Theta +
      1 / n * (term1 + term2 + term3 + term4)
  }

  Hess_Theta
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
  grad_Theta <- objfun(
    Ly, Ltid, Theta, lambda, sigma2, theta_mu, tau, "grad"
  )[['grad_Theta']]
  ThGZ <- as.matrix(t(Theta) %*% G %*% grad_Theta)
  DelGZ <- as.matrix(t(Delta) %*% G %*% grad_Theta)
  Dproj <- - 0.5 * Delta %*% (ThGZ + t(ThGZ)) -
    0.5 * Theta %*% (DelGZ + t(DelGZ))

  manifold.Stiefel.project(Hess_Theta, Theta, G) + Dproj
}

Delta <- matrix(rnorm(p*q), p, q)
Delta <- manifold.Stiefel.project(Delta, Theta, G)
Hess_Theta_apply(Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2)
maniHess_Theta_apply(Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2)

check_H_Theta <- check.derivatives(
  0,
  \(h) {
    grad_Theta <- objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta + h * Delta, lambda, sigma2
    )[["grad_Theta"]]
    grad_Theta
  },
  \(h) {
    Hess_Theta_apply(
      Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    ) |> as.vector()
  }
)

grad_Theta <- objfun(
  dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
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
    grad_Theta <- objfun(
      dat$Ly[1:N], dat$Ltid[1:N], Theta + h * Delta, lambda, sigma2
    )[["grad_Theta"]]
    manifold.Stiefel.project(grad_Theta, Theta + h * Delta, G) |> as.vector()
  },
  \(h) {
    maniHess_Theta_apply(
      Delta, dat$Ly[1:N], dat$Ltid[1:N], Theta, lambda, sigma2
    ) |> as.vector()
  }
)
