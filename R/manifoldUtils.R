sym <- function(Z, X, G = NULL) {
  p <- nrow(Z)
  q <- ncol(Z)
  stopifnot(p == nrow(X) && q == ncol(X))
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
  }
  return(sum(X * G %*% Y))
}

manifold.Stiefel.isIn <- function(X, G = NULL) {
  q <- ncol(X)
  if (is.null(G)) {
    return(all.equal(as.matrix(crossprod(X)), diag(1, q, q)))
  }
  all.equal(as.matrix(t(X) %*% G %*% X), diag(1, q, q))
}

manifold.Stiefel.isInTangent <- function(Z, X, G = NULL) {
  q <- ncol(X)
  if (is.null(G)) {
    return(all.equal(
      as.matrix(crossprod(X, Z)),
      -as.matrix(crossprod(Z, X))
    ))
  }
  all.equal(
    as.matrix(crossprod(X, G %*% Z)),
    as.matrix(-crossprod(Z, G %*% X))
  )
}

manifold.Stiefel.project <- function(Z, X, G = NULL) {
  # project Z onto the tangent space at X
  symZGX <- sym(Z, X, G)
  Z <- Z - X %*% symZGX
  return(as.matrix(Z))
}

manifold.Stiefel.transport <- function(Z, Y, X = NULL, G = NULL) {
  # transport Z from T_X to T_Y
  # TODO: find an isometric vector transport
  manifold.Stiefel.project(Z, Y, G)
}

manifold.Stiefel.retract <- function(Z, X = NULL, G = NULL) {
  # Retr_X(Z)
  # If X == NULL, view Z as X + Z.
  XZ <- if (is.null(X)) Z else X + Z
  if (is.null(G)) {
    M <- as.matrix(Matrix::crossprod(XZ))
  } else {
    M <- as.matrix(t(XZ) %*% G %*% XZ)
  }
  R <- Matrix::chol(M)
  as.matrix(XZ %*% solve(R))
}

# Enforce consistent QR branch: flip columns of X so that diag(Y^T G X) > 0
orient_Y_wrt_X <- function(X, Y, G = NULL) {
  if (is.null(G)) {
    C <- crossprod(X, Y) # C = X^T Y  (r x r)
  } else {
    C <- crossprod(X, G %*% Y) # C = X^T G Y  (r x r)
  }
  s <- sign(diag(C))
  s[s == 0] <- 1
  list(Yo = sweep(Y, 2, s, "*"), Co = sweep(C, 2, s, "*"), sgn = s)
}

# Solve for the upper-triangular U in: C U + U^T C^T = 2 I  (unique with diag(U)>0 near I)
solve_upper_inv_rectract <- function(C) {
  r <- ncol(C)
  stopifnot(nrow(C) == r)
  U <- matrix(0.0, r, r)
  for (j in 1:r) {
    M <- C[1:j, 1:j, drop = FALSE] # j x j
    rhs <- numeric(j)
    if (j > 1) {
      for (i in 1:(j - 1)) {
        # rhs_i = - sum_{k=1}^i C[j, k] * U[k, i]  (correct formula)
        rhs[i] <- -sum(C[j, 1:i] * U[1:i, i])
      }
    }
    rhs[j] <- 1.0
    u <- qr.solve(M, rhs) # robust small dense solve
    U[1:j, j] <- u
    # ensure strictly positive diagonals to stick to the same QR branch
    if (U[j, j] <= 0) U[j, j] <- abs(U[j, j]) + .Machine$double.eps
  }
  U
}

manifold.Stiefel.invRetract <- function(Y, X, G = NULL) {
  # Compute R_X^{-1}(Y)
  if (is.null(G)) {
    # X, Y: n×r with t(X) %*% X = I = t(Y) %*% Y
    oy <- orient_Y_wrt_X(X, Y) # G = I case
    U <- solve_upper_inv_rectract(oy$Co)
    Xi <- oy$Yo %*% U - X
    # attr(Xi, "tangent_violation_F") <- norm(crossprod(X, Xi) + crossprod(Xi, X), "F")
    return(Xi)
  }
  # Either supply G (spd) or its Cholesky R with G = t(R) %*% R
  oy <- orient_Y_wrt_X(X, Y, G = G) # Co = X^T G Y with oriented columns
  U <- solve_upper_inv_rectract(oy$Co)
  Xi <- oy$Yo %*% U - X # tangent at X
  # tangent check: X^T G Xi + Xi^T G X
  # GXi <- G %*% Xi          # G %*% Xi without forming G
  # viol <- norm(crossprod(X, GXi) + crossprod(GXi, X), "F")
  # attr(Xi, "tangent_violation_F") <- viol
  Xi
}

manifold.Stiefel.normalize <- function(x, X0 = NULL, G = NULL) {
  # A intuitive normalization algorithm
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
