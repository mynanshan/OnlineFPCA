### computational helpers

### eigen and inverse --------------------

decomp_Mlambda <- function(K, R, lambda = 0, Rtype = c("diag", "general")) {
  ## Decompose the M(lambda) for penalized regression
  # with multiple smoothing parameters
  # K: a p.s.d. symmetric matrix, cov of the design
  # R: 1) length p vector, represents a diagonal penalty matrix
  #    2) general p.s.d. matrix
  # M(lambda) = K + lambda * R
  Rtype <- match.arg(Rtype)
  p <- nrow(K)
  stopifnot(ncol(K) == p)
  nlam <- length(lambda)
  if (Rtype == "diag") {
    R <- as.vector(R)
    stopifnot(length(R) == p)
    stopifnot(all(R > 0))
    Rsqrt <- sqrt(R)
    A <- sweep(sweep(K, 1, Rsqrt^(-1), "*"), 2, Rsqrt^(-1), "*")
    eig_obj <- eigen(A, symmetric = TRUE) # A =  U D U^T
    U <- eig_obj$vectors
    D <- eig_obj$values
    D_plus_lambda <- outer(D, lambda, "+")
    V <- sweep(U, 1, Rsqrt, "*")
    Vinv <- sweep(U, 1, Rsqrt^(-1), "*")
    return(list(V = V, Vinv = Vinv, diag_plus_lambda = D_plus_lambda))
  }
}

# X <- matrix(runif(70), 10, 7)
# Mdecomp <- decomp_Mlambda(t(X) %*% X, 1:7, 10^((-1):(-4)))
# (M1 <- (t(X) %*% X + 0.01 * diag(1:7)))
# Mdecomp$V %*% diag(Mdecomp$diag_plus_lambda[,2]) %*% t(Mdecomp$V)
# solve(M1)
# Mdecomp$Vinv %*% diag((Mdecomp$diag_plus_lambda[,2])^(-1)) %*% t(Mdecomp$Vinv)

eigen_lrpsd <- function(M) {
  ## eigendecomposition of p.s.d matrices with small eigenvalues removed
  eig_obj <- eigen(M, symmetric = TRUE)
  V <- eig_obj$vectors
  D <- eig_obj$values
  eigval_tol <- max(D) * .Machine$double.eps
  rm_id <- which(D < eigval_tol)
  if (length(rm_id) > 0) {
    V <- V[, -rm_id, drop = FALSE]
    D <- D[-rm_id]
  }
  return(list(vectors = V, values = D))
}


### array operations --------------------

asl <- function(a, i) {
  # asl for array slicing
  # 3d array only
  dims <- dim(a)
  return(array(a[,, i], dim = dims[-3]))
}


### grid points -----------------

margins2grid <- function(mlist) {
  ndim <- length(mlist)
  m <- sapply(mlist, length)
  grid <- matrix(nrow = prod(m), ncol = ndim)
  for (j in 1:ndim) {
    nrep_each <- ifelse(j > 1, prod(m[1:(j - 1)]), 1)
    nrep_whole <- ifelse(j < ndim, prod(m[(j + 1):ndim]), 1)
    grid[, j] <- mlist[[j]] |>
      (\(x) rep(x, each = nrep_each))() |>
      (\(x) rep(x, nrep_whole))()
  }
  return(grid)
}

### manipulating covariances --------------------

cov.form_covdat <- function(Lt, Ly = NULL, diag.rm = TRUE) {
  ## compute raw covariances
  ## form longitudinal data to a dataframe of design and responses
  if (is.list(Ly) && is.list(Lt)) {
    n <- length(Ly)
    stopifnot(length(Lt) == n)
    Lmi <- sapply(Lt, nrow)
    D <- ncol(Lt[[1]])
    covdat <- lapply(1:n, \(i) {
      stopifnot(nrow(Lt[[i]]) == length(Ly[[i]]))
      dat <- matrix(NA, Lmi[i]^2, 2 * D + 1)
      dat[, 2 * D + 1] <- rep(Ly[[i]], times = Lmi[i]) *
        rep(Ly[[i]], each = Lmi[i])
      dat[, 1:D] <- apply(Lt[[i]], 2, \(x) rep(x, times = Lmi[i]))
      dat[, (D + 1):(2 * D)] <- apply(Lt[[i]], 2, \(x) rep(x, each = Lmi[i]))
      if (diag.rm) {
        rm_id <- (0:(Lmi[i] - 1)) * Lmi[i] + 1:Lmi[i]
        dat <- dat[-rm_id, ]
      }
      return(dat)
    })
    covdat <- data.frame(do.call(rbind, covdat))
    names(covdat) <- c(paste0("s", 1:D), paste0("t", 1:D), "z")
  } else if (is.matrix(Lt)) {
    m <- nrow(Lt)
    D <- ncol(Lt)
    covdat <- matrix(NA, m^2, 2 * D)
    covdat[, 1:D] <- apply(Lt, 2, \(x) rep(x, times = m))
    covdat[, (D + 1):(2 * D)] <- apply(Lt, 2, \(x) rep(x, each = m))
    covdat <- data.frame(covdat)
    names(covdat) <- c(paste0("s", 1:D), paste0("t", 1:D))
    if (!is.null(Ly) && is.numeric(Ly)) {
      stopifnot(length(Ly) == m)
      covdat <- cbind(covdat, z = rep(Ly, times = m) * rep(Ly, each = m))
    }
    if (diag.rm) {
      rm_id <- (0:(m - 1)) * m + 1:m
      covdat <- covdat[-rm_id, ]
    }
  } else {
    stop("Wrong input types!")
  }
  return(covdat)
}

cov.form_diagdat <- function(Lt, Ly = NULL) {
  if (is.list(Ly) && is.list(Lt)) {
    stopifnot(length(Lt) == length(Ly))
    D <- ncol(Lt[[1]])
    dat <- data.frame(cbind(do.call(rbind, Lt), unlist(Ly)^2))
    names(dat) <- c(paste0("t", 1:D), "z")
  } else if (is.matrix(Lt)) {
    D <- ncol(Lt)
    dat <- data.frame(cbind(Lt, Ly^2))
    names(dat) <- paste0("t", 1:D)
    if (!is.null(Ly) && is.numeric(Ly)) {
      dat <- cbind(dat, z = Ly^2)
    }
  }
  return(dat)
}


cov.symmetrize <- function(C) {
  stopifnot(is.matrix(C))
  stopifnot(nrow(C) == ncol(C))
  C <- (C + t(C)) / 2
  return(C)
}


## conversion between longitudinal data and dataframes --------
df2Ldata = function(Tmat, Yvec, Lmi) {
  N = length(Lmi)
  Tmat = as.matrix(Tmat)
  subj_start_id = c(1, cumsum(Lmi[-N]) + 1)
  subj_end_id = cumsum(Lmi)
  Ly = lapply(
    seq_len(N),
    \(i) {
      Yvec[subj_start_id[i]:subj_end_id[i]]
    }
  )
  Lt = lapply(
    seq_len(N),
    \(i) {
      Tmat[subj_start_id[i]:subj_end_id[i], ]
    }
  )
  return(list(Ly = Ly, Lt = Lt))
}


## manipulating eigenfunctions ----------
flip_direc <- function(Phi_target, Phi_ref) {
  idx <- (colMeans(Phi_target * Phi_ref) < 0)
  Phi_target[, idx] <- -Phi_target[, idx]
  attr(Phi_target, "flipped") <- idx
  return(Phi_target)
}

match_fpc <- function(Phi_target, Phi_ref) {
  npc <- ncol(Phi_target)
  stopifnot(ncol(Phi_ref) == npc)
  if (npc > 1) {
    match_id <- numeric(npc)
    for (k in 1:npc) {
      err <- pmin(
        colMeans(sweep(Phi_target, 1, Phi_ref[, k], "-")^2),
        colMeans(sweep(Phi_target, 1, Phi_ref[, k], "+")^2)
      )
      ord <- order(err, decreasing = FALSE)
      for (j in 1:npc) {
        if (ord[j] %in% match_id) {
          next
        } else {
          match_id[k] <- ord[j]
          if (err[match_id[k]] > 0.1 * min(err[-match_id[k]])) {
            warning("Probably a bad fit.")
          }
          break
        }
      }
    }
    if (any(table(match_id) > 1)) {
      stop("Wrong match.")
    }
    Phi_target <- Phi_target[, match_id]
    attr(Phi_target, "match_id") <- match_id
  } else {
    attr(Phi_target, "match_id") <- 1
  }
  return(flip_direc(Phi_target, Phi_ref))
}

## Bridges to {funData} -------------------
# gridData2funData <- function(tlist, Y) {
#   d <- sapply(tlist, length)
#   N <- nrow(Y)
#   p <- ncol(Y)
#   stopifnot(p == prod(d))
#   return(funData::funData(argvals = tlist, X = array(Y, dim = c(N, d))))
# }
