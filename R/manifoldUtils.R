sym <- function(Z, X, G = NULL) {
  p <- nrow(Z); q <- ncol(Z)
  stopifnot(p == nrow(X) && q ==ncol(X))
  if (is.null(G)) {
    Zt_X <- Matrix::crossprod(Z, X)
  } else {
    Zt_X <- Matrix::crossprod(Z, G %*% X)
  }
  0.5 * as.matrix(Zt_X + Matrix::t(Zt_X))
}

manifold.Stiefel.inprod <- function(X, Y, G = NULL) {
  if (is.null(G)) {
    return(sum(X * Y))
  } else {
    return(sum(X * G %*% Y))
  }
}

manifold.Stiefel.project <- function(Z, X, G = NULL) {
  # project Z onto the tangent space at X
  symZGX <- sym(Z, X, G)
  Z <- Z - X %*% symZGX
  return(as.matrix(Z))
}

manifold.Stiefel.retract <- function(X, G = NULL) {
  if (is.null(G)) {
    M <- as.matrix(Matrix::crossprod(X))
  } else {
    M <- as.matrix(t(X) %*% G %*% X)
  }
  R <- Matrix::chol(M)
  X <- as.matrix(X %*% solve(R))
  return(X)
}

manifold.Stiefel.retractSinglePC <- function(x, X0 = NULL, G = NULL) {
  # Suppose X = [X0, x]
  # X0 part satisfies the manifold constraint, but x part does not
  # Now project x s.t. X is retracted onto the Stiefel manifold
  x <- as.matrix(x)
  # stopifnot(ncol(x)==1)
  if (!is.null(X0)) {
    X0 <- as.matrix(X0)
    # stopifnot(nrow(x)==nrow(X0))
    if (is.null(G)) {
      inprod_X0_x <- crossprod(X0, x)
    } else {
      # stopifnot(nrow(x)==nrow(X0))
      inprod_X0_x <- as.matrix(crossprod(X0, G %*% x))
    }
    x <- x - X0 %*% inprod_X0_x
  }
  x <- as.vector(x)
  if (is.null(G)) {
    x <- x / sqrt(sum(x^2))
  } else {
    x <- x / sqrt(sum(x * G %*% x))
  }
  return(x)
}