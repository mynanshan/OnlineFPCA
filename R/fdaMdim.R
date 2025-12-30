# require(fda)
# require(rTensor)
# require(mgcv)

library(Matrix)

### Extensions to {fda} =================

## correspondence between basis2d and coef:
## basisList = list(basis1, basis2, ..., basisD)
## dim(coef) = c(lD, ..., l2, l1)

TensorBasis <- function(basislist) {
  # tensor product basis
  D <- length(basislist)
  nbasis <- 1
  for (d in 1:D) {
    stopifnot(is(basislist[[d]], "basisfd"))
    nbasis <- nbasis * basislist[[d]]$nbasis
  }
  basis <- structure(basislist, class = "TensorBasis")
  attr(basis, "dimension") <- D
  attr(basis, "nbasis") <- nbasis
  return(basis)
}


FuncData <- function(coef, basis) {
  # multidimensional functional data
  if (!is.matrix(coef)) {
    coef <- as.matrix(coef)
  }
  if (is(basis, "basisfd")) {
    return(fda::fd(coef, basis))
  } else if (is(basis, "TensorBasis")) {
    nb <- sapply(basis, \(x) x[["nbasis"]])
    ntb <- prod(nb)
    stopifnot(nrow(coef) == ntb)
    fdobj <- structure(list(coefs = coef, basis = basis), class = "FuncData")
    attr(fdobj, "dimension") <- attr(basis, "dimension")
    return(fdobj)
  } else {
    stop("class(basis) need to be basisfd or tensorBasis")
  }
}


eval_basis <- function(evalarg, basis, Lfdobj = 0, sparse_matrix = TRUE) {
  if (is(basis, "basisfd")) {
    evalarg <- as.vector(evalarg)
    if (length(Lfdobj) > 1) {
      Lfdobj <- Lfdobj[1]
      warning("Only the 1st element of Lfdobj is used.")
    }
    basismat <- fda::eval.basis(evalarg, basis, Lfdobj)
    if (sparse_matrix && basis$type == "bspline") {
      basismat <- as(basismat, "sparseMatrix")
    }
    return(basismat)
  } else if (is(basis, "TensorBasis")) {
    evalarg <- as.matrix(evalarg)
    nt <- nrow(evalarg)
    D <- length(basis)
    stopifnot(ncol(evalarg) == D)
    if (length(Lfdobj) == 1) {
      Lfdobj <- rep(Lfdobj, D)
    }
    basismat_margin <- lapply(
      1:D,
      \(d) {
        B <- fda::eval.basis(evalarg[, d], basis[[d]], Lfdobj[d])
        if (sparse_matrix && basis[[d]]$type == "bspline") {
          B <- as(B, "sparseMatrix")
        }
        return(B)
      }
    )
    basismat <- mgcv::tensor.prod.model.matrix(basismat_margin)
    return(basismat)
  } else {
    stop("class(basis) need to be basisfd or tensorBasis")
  }
}


eval_fd <- function(evalarg, fdobj, Lfdobj = 0) {
  if (is(fdobj, "fd")) {
    evalarg <- as.vector(evalarg)
    if (length(Lfdobj) > 1) {
      Lfdobj <- Lfdobj[1]
      warning("Only the 1st element of Lfdobj is used.")
    }
    return(fda::eval.fd(evalarg, fdobj, Lfdobj))
  } else if (is(fdobj, "FuncData")) {
    evalarg <- as.matrix(evalarg)
    nt <- nrow(evalarg)
    D <- length(fdobj$basis)
    if (length(Lfdobj) == 1) {
      Lfdobj <- rep(Lfdobj, D)
    }
    basismat <- eval_basis(evalarg, fdobj$basis, Lfdobj)
    return(as.matrix(basismat %*% fdobj$coef))
  } else {
    stop("class(fdobj) need to be fd or FuncData")
  }
}


get_basis_inprod_matrix <- function(basis, sparse_matrix = TRUE) {
  if (is(basis, "basisfd")) {
    inprodmat <- fda::inprod(basis, basis)
    if (sparse_matrix) {
      inprodmat <- as(inprodmat, "sparseMatrix")
    }
    return(inprodmat)
  } else if (is(basis, "TensorBasis")) {
    D <- length(basis)
    inprodmat_margin <- lapply(
      1:D,
      \(d) {
        G <- fda::inprod(basis[[d]], basis[[d]])
        if (sparse_matrix && basis[[d]]$type == "bspline") {
          G <- as(G, "sparseMatrix")
        }
        return(G)
      }
    )
    inprodmat <- inprodmat_margin[[1]]
    if (D > 1) {
      for (d in 2:D) {
        inprodmat <- inprodmat %x% inprodmat_margin[[d]]
      }
    }
    return(inprodmat)
  } else {
    stop("class(basis) need to be basisfd or tensorBasis")
  }
}


get_basis_penalty_matrix <- function(
  basis,
  penLfd = 2,
  sparse_matrix = TRUE,
  sum_matrices = TRUE
) {
  if (is(basis, "basisfd")) {
    penmat <- fda::inprod(basis, basis, penLfd, penLfd)
    if (sparse_matrix && basis$type == "bspline") {
      penmat <- as(penmat, "sparseMatrix")
    }
    return(penmat)
  } else if (is(basis, "TensorBasis")) {
    D <- length(basis)
    penmat_margin <- lapply(
      1:D,
      \(d) {
        S <- fda::inprod(basis[[d]], basis[[d]], penLfd, penLfd)
        if (sparse_matrix && basis[[d]]$type == "bspline") {
          S <- as(S, "sparseMatrix")
        }
        return(S)
      }
    )
    penmat_list <- mgcv::tensor.prod.penalties(penmat_margin)
    if (sum_matrices) {
      penmat <- penmat_list[[1]]
      if (D > 1) {
        for (d in 2:D) {
          penmat <- penmat + penmat_list[[d]]
        }
      }
      return(penmat)
    }
  } else {
    stop("class(basis) need to be basisfd or tensorBasis")
  }
}


smooth_basis <- function(argvals, Y, basis, lambda = 0, penLfd = 2) {
  B <- eval_basis(argvals, basis)
  S <- get_basis_penalty_matrix(basis, penLfd, sum_matrices = TRUE)
  coef <- as.matrix(Matrix::solve(
    Matrix::crossprod(B) + lambda * S,
    Matrix::crossprod(B, Y)
  ))
  if (is(basis, "basisfd")) {
    return(fda::fd(coef = coef, basisobj = basis))
  } else if (is(basis, "TensorBasis")) {
    return(FuncData(coef, basis))
  } else {
    stop("class(basis) need to be basisfd or tensorBasis")
  }
}


smooth_basis.gcv <- function(argvals, Y, basis, lambda = 0, penLfd = 2) {
  B <- eval_basis(argvals, basis) # = B1 %x% B2
  S <- get_basis_penalty_matrix(basis, penLfd, sum_matrices = TRUE)
  BtB <- Matrix::crossprod(B)
  BtY <- Matrix::crossprod(B, Y)

  # Bmat_list <- lapply(basis, \(b) eval_basis(argvals, basis))
  # Smat_list <- lapply(basis, \(b) fda::inprod(basis, basis, penLfd, penLfd))
  # eigenB_list <- lapply(Bmat_list, \(Bi)
  #                       eigen(Matrix::crossprod(Bi), symmetric = TRUE))

  nt <- nrow(B)
  nlambda <- length(lambda)
  coef <- matrix(nrow = ncol(B), ncol = nlambda)
  MSE <- numeric(nlambda)
  den <- numeric(nlambda)
  GCV <- numeric(nlambda)
  for (i in 1:nlambda) {
    Minv <- Matrix::solve(BtB + lambda[i] * S)
    coef[, i] <- as.matrix(Minv %*% BtY)
    Ypred <- B %*% coef[, i]
    MSE[i] <- mean((Y - Ypred)^2)
    # smoothing matrix: S = Psi %*% Minv %*% t(Psi)
    den[i] <- mean(1 - sapply(1:nt, \(j) sum(B[j, ] * Minv %*% B[j, ])))^2
    GCV[i] <- MSE[i] / den[i]
  }
  idx <- which.min(GCV)
  if (is(basis, "basisfd")) {
    fdObj <- fda::fd(coef = coef[, idx], basisobj = basis)
  } else if (is(basis, "TensorBasis")) {
    fdObj <- FuncData(coef[, idx], basis)
  } else {
    stop("class(basis) need to be basisfd or tensorBasis")
  }
  return(list(
    lambda = lambda,
    GCV = GCV,
    coef = coef,
    lambda.min = lambda[idx],
    fd = fdObj
  ))
}


bind_argvals <- function(evalArgs, is.basisfd = TRUE) {
  # evalArgs: a list of evaluation arguments
  if (is.basisfd) {
    return(do.call(c, evalArgs))
  } else {
    return(do.call(rbind, evalArgs))
  }
}
