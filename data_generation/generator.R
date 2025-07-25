### functional data generator

library(fda)


# Random process settings ============================

get_rp.low1 <- function() {
  t0 <- 0; t1 <- 1
  mu <- function(x) rep(0, length(x))
  eigfun <- function(x) {
    return(cbind(-cos(pi*x)*sqrt(2), sin(pi*x)*sqrt(2)))
  }
  eigval <- c(4,1)
  return(RandomProcess(
    ndim = 1, t0 = t0, t1 = t1, meanfun = mu,
    eigfun = eigfun, eigval = eigval, type = "gaussian", rank = 2
  ))
}

get_rp.xiao2018b <- function(nu=1, phi=0.07, scale=1) {
  t0 <- 0; t1 <- 1
  mu <- function(x) rep(0, length(x))
  covfun <- function(x) {
    distmat <- as.matrix(dist(x,diag=T,upper=T))
    return(scale * kernel.matern(distmat,nu,phi))
  }
  return(RandomProcess(
    ndim = 1, t0 = t0, t1 = t1, meanfun = mu,
    covfun = covfun, type = "gaussian", rank = Inf
  ))
}

get_rp.yao2005 <- function(demean = FALSE) {
  t0 <- 0; t1 <- 10
  if (demean) {
    mu <- function(x) rep(0, length(x))
  } else {
    mu <- function(x) x + sin(x)
  }
  eigfun <- function(x) {
    return(cbind(-cos(pi*x/10)/sqrt(5), sin(pi*x/10)/sqrt(5)))
  }
  eigval <- c(4,1)
  return(RandomProcess(
    ndim = 1, t0 = t0, t1 = t1, meanfun = mu,
    eigfun = eigfun, eigval = eigval, type = "gaussian", rank = 2
  ))
}

get_rp.yang2021 <- function(alpha = 2, npc = 10) {
  t0 <- 0; t1 <- 1
  stopifnot(npc > 1)
  # mu <- function(x) 2 * sin(2*pi*x)
  mu <- function(x) rep(0, length(x))
  lambda0 <- 0.4
  eigval <- lambda0 * (1:npc)^(-alpha)
  eigsum <- sum(0.4 * (1:10)^(-2))
  eigval <- eigval / sum(eigval) * eigsum
  eigfun <- function(x) {
    x <- as.numeric(x)
    n <- length(x)
    efun1 <- rep(1, n)
    efun <- matrix(sqrt(2) * sapply(2:npc, \(i) cos((i-1)*pi*x)), n, npc-1)
    return(unname(cbind(efun1, efun)))
  }
  return(RandomProcess(
    ndim = 1, t0 = t0, t1 = t1, meanfun = mu,
    eigfun = eigfun, eigval = eigval, type = "gaussian", rank = npc
  ))
}

get_rp.he2022 <- function(set = c("easy", "prac")) {
  set <- match.arg(set)
  t0 <- 0; t1 <- 1
  mu <- function(x) rep(0, length(x))
  if (set == "easy") {
    R <- 5
    npc <- 3
    eigval <- c(1, 0.66, 0.517)
  }
  if (set == "prac") {
    R <- 10
    npc <- 5
    eigval <- c(1, 0.66, 0.517, 0.435, 0.381)
  }
  Q <- pracma::randortho(R)[,1:npc]
  eigfun <- local(function(x) {
    x <- as.numeric(x)
    n <- length(x)
    efun <- matrix(sqrt(2) * sapply(1:R, \(i) sin(i*pi*x)), n, R)
    return(efun %*% Q)
  })
  return(RandomProcess(
    ndim = 1, t0 = t0, t1 = t1, meanfun = mu,
    eigfun = eigfun, eigval = eigval, type = "gaussian", rank = npc
  ))
}


get_rp.wang2020 <- function(r1 = 3, r2 = 2, alpha = 2, lambda0 = 1) {
  npc <- r1 * r2
  t0 <- c(0, 0); t1 <- c(1, 1)
  mu <- function(X) {return(rep(0, nrow(X)))}
  eigval <- lambda0 * exp(-(0:(npc-1)) * log(alpha))
  eigfun <- local({
    efun <- function(x, r) {
      return(matrix(sapply(1:r, \(k) sqrt(2)*cos(pi*x*k)), length(x), r))
    }
    function(X) {
      N <- nrow(X)
      psi1 <- efun(X[,1], r1)
      psi2 <- efun(X[,2], r2)
      phimat <- matrix(nrow = N, ncol = npc)
      for (k in 1:npc) {
        i <- (k-1) %/% r2 + 1
        j <- (k-1) %% r2 + 1
        phimat[,k] <- psi1[,i] * psi2[,j]
      }
      return(phimat)
    }
  })
  return(RandomProcess(
    ndim = 2, t0 = t0, t1 = t1, meanfun = mu,
    eigfun = eigfun, eigval = eigval, type = "gaussian", rank = npc
  ))
}


RandomProcess <- function(meanfun, covfun = NULL,
                          eigfun = NULL, eigval = NULL,
                          ndim = 1, t0, t1, rank = Inf,
                          type = c("gaussian")) {
  type <- match.arg(type)
  out <- list(
    t0 = t0, t1 = t1, ndim = ndim, meanfun = meanfun,
    covfun = covfun, eigfun = eigfun, eigval = eigval,
    type = type, rank = rank)
  return(out)
}



# Functions for data generation ============================

get_measurements <-
  function(rp, n,
           m = NULL, m_min = NULL, m_max = NULL, m_mean = NULL, m_sd = NULL,
           neval = 51, design_type = c("fixed", "random"),
           m_type = c("unif", "gaussian"), sigma = 0.5) {
  design_type <- match.arg(design_type)
  m_type <- match.arg(m_type)
  if (design_type == "random") {
    if (m_type == "unif") {
      if (!is.null(m) && (is.null(m_min) || is.null(m_max))) {
        m_min <- m_max <- m[1]
      }
      stopifnot(!is.null(m_min) && !is.null(m_max))
      stopifnot(m_max >= m_min)
    }
    if (m_type == "gaussian") {
      if (!is.null(m) && is.null(m_mean)) {
        m_mean <- m[1]
      }
      if (is.null(m_sd)) m_sd <- 1
      stopifnot(!is.null(m_mean) && !is.null(m_sd))
    }
  } else {
    stopifnot(!is.null(m))
    if (rp$ndim > 1 && length(m) == 1) {
      m <- rep(m, rp$ndim)
    }
    stopifnot(length(m) == rp$ndim)
  }
  lowrank <- !is.infinite(rp$rank)
  if (lowrank) {
    q <- rp$rank
    if (rp$type == "gaussian") {
      scores <- matrix(rnorm(n * q, mean = 0, sd = 1), n, q)
      scores <- sweep(scores, 2, sqrt(rp$eigval), "*")
    } else {
      stop("NOT IMPLEMENTED")
    }
  }
  
  if (design_type == "fixed") {
    if (rp$ndim == 1) {
      tgrid_margin <- list(seq(rp$t0, rp$t1, length.out = m))
      tgrid <- matrix(tgrid_margin[[1]], ncol = 1)
    } else {
      tgrid_margin <- lapply(1:rp$ndim, \(j) {
        seq(rp$t0[j], rp$t1[j], length.out = m[j])
      })
      tgrid <- matrix(nrow = prod(m), ncol = rp$ndim)
      for (j in 1:rp$ndim) {
        # nrep_each <- ifelse(j < rp$ndim, prod(m[(j+1):rp$ndim]), 1)
        # nrep_whole <- ifelse(j > 1, prod(m[1:(j-1)]), 1)
        nrep_each <- ifelse(j > 1, prod(m[1:(j-1)]), 1)
        nrep_whole <- ifelse(j < rp$ndim, prod(m[(j+1):rp$ndim]), 1)
        tgrid[,j] <- tgrid_margin[[j]] |>
          (\(x) rep(x, each = nrep_each))() |>
          (\(x) rep(x, nrep_whole))()
      }
    }
    ntgrid <- prod(m)
    mu_t <- rp$meanfun(tgrid)
    if (lowrank) { # generate by eigenfunctions
      Phi_t <- rp$eigfun(tgrid)
      X = Phi_t %*% t(scores)
      Y <- sweep(Phi_t %*% t(scores), 1, mu_t, "+") +
        rnorm(n * ntgrid, mean = 0, sd = sigma)
      return(list(tgrid = tgrid, X = X, Y = Y, scores = scores, tgrid_margin = tgrid_margin))
    } else { # generate by covariance
      covmat <- rp$covfun(tgrid)
      if (rp$type == "gaussian") {
        Z <- matrix(rnorm(n * ntgrid), nrow = n, ncol = ntgrid)
      } else {stop("NOT IMPLEMENTED")}
      diag(covmat) <- diag(covmat) + 1e-12
      A <- chol(covmat)
      Y <- sweep(Z %*% A, 2, mu_t, "+") +
        rnorm(n * ntgrid, mean = 0, sd = sigma)
      return(list(tgrid = tgrid, Y = Y, tgrid_margin = tgrid_margin))
    }
  } else if (design_type == "random") {
    if (length(neval) == 1 && rp$ndim > 1) neval <- rep(neval, rp$ndim)
    stopifnot(length(neval) == rp$ndim)
    if (rp$ndim == 1) {
      tgrid <- seq(rp$t0, rp$t1, length.out = neval)
    } else {
      tgrid_margin <- lapply(1:rp$ndim, \(j) {
        seq(rp$t0[j], rp$t1[j], length.out = neval[j])
      })
      tgrid <- matrix(nrow = prod(neval), ncol = rp$ndim)
      for (j in 1:rp$ndim) {
        # nrep_each <- ifelse(j < rp$ndim, prod(neval[(j+1):rp$ndim]), 1)
        # nrep_whole <- ifelse(j > 1, prod(neval[1:(j-1)]), 1)
        nrep_each <- ifelse(j > 1, prod(neval[1:(j-1)]), 1)
        nrep_whole <- ifelse(j < rp$ndim, prod(neval[(j+1):rp$ndim]), 1)
        tgrid[,j] <- tgrid_margin[[j]] |>
          (\(x) rep(x, each = nrep_each))() |>
          (\(x) rep(x, nrep_whole))()
      }
    }
    ntgrid <- prod(neval)
    Ly <- vector("list", n)
    Lt <- vector("list", n)
    Ltid <- vector("list", n)
    if (m_type == "unif") {
      r <- m_max - m_min + 1
      Lmi <- floor(runif(n) * r) + m_min
    }
    if (m_type == "gaussian") {
      Lmi <- pmax(round(rnorm(n, mean = m_mean, sd = m_sd)), 1)
    }
    Lmi = pmin(Lmi, neval)
    for (i in 1:n) {
      # print(Lmi[i])
      Ltid[[i]] <- idx <- sort(sample(
        1:ntgrid, Lmi[i], replace = FALSE))
      if (rp$ndim == 1) {
        Lt[[i]] <- Ti <- tgrid[idx]
      } else {
        Lt[[i]] <- Ti <- tgrid[idx,,drop=F]
      }
      mu_ti <- rp$meanfun(Ti)
      if (lowrank) { # generate by eigenfunctions
        Phi_ti <- rp$eigfun(Ti)
        Ly[[i]] <- as.vector(mu_ti + Phi_ti %*% scores[i,]) +
          rnorm(Lmi[i], mean = 0, sd = sigma)
      } else { # generate by covariance
        covmat <- rp$covfun(Ti)
        diag(covmat) <- diag(covmat) + 1e-12
        if (rp$type == "gaussian") {
          z <- rnorm(Lmi[i])
        } else {stop("NOT IMPLEMENTED")}
        A <- chol(covmat)
        Ly[[i]] <- as.vector(
          t(A) %*% z + mu_ti + rnorm(Lmi[i], mean = 0, sd = sigma))
      }
    }
    out <- list(Lt = Lt, Ly = Ly, Lmi = Lmi, Ltid = Ltid, tgrid = tgrid)
    if (lowrank) out[["scores"]] <- scores
    return(out)
  }
}


# Measurements from real eigenfunctions  ----------------------------------

gendata.gfr <- function(
  n = 1000, t0 = 0, t1 = 1,
  m_mean = 6, m_sd = 2, sigma = 0.08,
  eigf.path = "application/eigf_gfr.rds",
  eigval.path = "application/eigval_gfr.rds"
) {
  PhiMat <- readRDS(eigf.path)  # (M, npc)
  lamvec <- readRDS(eigval.path)
  npc <- ncol(PhiMat)
  stopifnot(length(lamvec) == npc)
  M <- nrow(PhiMat)
  tgrid <- seq(t0, t1, length.out = M)
  scores <- rnorm(n * npc) |> 
    matrix(nrow = n, ncol = npc) |> 
    sweep(2, sqrt(lamvec), "*")
  X <- scores %*% t(PhiMat)  # (n, M)
  Y <- X + rnorm(n * M, mean = 0, sd = sigma)
  Lmi <- rnorm(n, mean = m_mean, sd = m_sd) |>
    pmax(2) |> round()
  Ltid <- lapply(1:n, \(i) sort(sample(1:M, Lmi[i])))
  Lt <- lapply(1:n, \(i) tgrid[Ltid[[i]]])
  Ly <- lapply(1:n, \(i) Y[i, Ltid[[i]]])
  list(
    Lt = Lt, Ly = Ly, Lmi = Lmi,
    Ltid = Ltid, tgrid = tgrid, t0 = t0, t1 = t1,
    Phi = PhiMat, lam = lamvec, sigma2 = sigma^2
  )
}

gendata.aqi <- function(
  n = 1000, t0 = c(0, 0), t1 = c(1, 1),
  m_mean = 30, m_sd = 6, sigma = 0.1,
  eigf.path = "application/eigf_aqi.rds",
  eigval.path = "application/eigval_aqi.rds",
  tgrid.path = "application/tgrid_aqi.rds"
) {
  npc <- 3
  PhiMat <- readRDS(eigf.path)[,1:npc]  # (M, npc)
  # PhiMat <- readRDS(eigf.path)[,c(1,2,4)]  # (M, npc)
  lamvec <- readRDS(eigval.path)[1:npc]  # (M, npc)
  # lamvec <- readRDS(eigval.path)[c(1,2,4)]  # (M, npc)
  # lamvec <- lamvec / mean(lamvec)
  tgrid <- readRDS(tgrid.path)
  # alpha <- 2
  # lamvec <- 4 * exp(-(0:(npc-1)) * log(alpha))
  # lamvec <- 4 * (1:npc)^(-alpha) 
  M <- nrow(PhiMat)
  stopifnot(length(lamvec) == npc)
  stopifnot(M == nrow(tgrid))
  stopifnot(ncol(tgrid) == 2)
  scores <- rnorm(n * npc) |> 
    matrix(nrow = n, ncol = npc) |> 
    sweep(2, sqrt(lamvec), "*")
  X <- scores %*% t(PhiMat)  # (n, M)
  Y <- X + rnorm(n * M, mean = 0, sd = sigma)
  Lmi <- rnorm(n, mean = m_mean, sd = m_sd) |>
    pmax(2) |> round()
  Ltid <- lapply(1:n, \(i) sort(sample(1:M, Lmi[i])))
  Lt <- lapply(1:n, \(i) tgrid[Ltid[[i]],,drop=F])
  Ly <- lapply(1:n, \(i) Y[i, Ltid[[i]]])
  list(
    Lt = Lt, Ly = Ly, Lmi = Lmi,
    Ltid = Ltid, tgrid = tgrid, t0 = t0, t1 = t1,
    Phi = PhiMat, lam = lamvec, sigma2 = sigma^2
  )
}
