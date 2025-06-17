pspc_stats <- function(
    Ly, Lt, basis, Theta, theta_mu, sigma2, lambda,
    stats = c("grad", "mse", "mse_offdiag")) {
  # consts: storing constants p, q, n, m2, Lmi
  if (is(basis, "basisfd")) {
    p <- basis$nbasis
  } else if (is(basis, "TensorBasis")) {
    p <- attr(basis, "nbasis")
  }
  n <- length(Ly)
  stopifnot(length(Lt) == n)
  Lmi <- sapply(Ly, length)
  
  if (length(dim(Theta)) == 3) {
    multiSmoothness <- TRUE
    q <- dim(Theta)[2]
    ntau <- dim(Theta)[3]
  } else if (length(dim(Theta)) == 2) {
    multiSmoothness <- FALSE
    q <- ncol(Theta)
    ntau <- 1
    Theta <- array(Theta, dim = c(dim(Theta),1))
    lambda <- matrix(lambda, nrow = length(lambda), ncol = 1)
  }
  stopifnot(dim(Theta)[1] == p)
  stopifnot(dim(Theta)[2] == q)
  stopifnot(nrow(lambda) == q)
  stopifnot(length(sigma2) == ntau)
  stopifnot(ncol(lambda) == ntau)
  
  if ("grad" %in% stats) {
    grad_Theta <- array(0, dim = c(p, q, ntau))
    grad_zeta <- numeric(ntau)
    grad_eta <- array(0, dim = c(q, ntau))
  }
  if ("mse" %in% stats) {
    mse <- numeric(ntau)
  }
  if ("mse_offdiag" %in% stats) {
    mse_offdiag <- numeric(ntau)
  }
  
  out <- list()
  
  for (i in 1:n) {
    Yi <- Ly[[i]]; Ti <- Lt[[i]]; mi <- Lmi[i]
    Bi <- eval_basis(Ti, basis)
    mu_i <- as.matrix(Bi %*% theta_mu)
    Yibar <- Yi - mu_i
    Bt_Yi <- as.matrix(Matrix::crossprod(Bi, Yibar))  # p x 1
    Yit_Yi <- sum(Yibar^2)
    
    Phi_i <- array(apply(Theta, 3, \(x) as.matrix(Bi %*% x)), dim = c(mi,q,ntau))  # mi x q x ntau
    Phit_Phi <- array(apply(Phi_i, 3, \(X) crossprod(X)), dim = c(q,q,ntau))  # q x q x ntau
    Phit_Yi <- apply(Phi_i, 3, \(X) crossprod(X, Yibar))  # q x ntau
    
    if ("grad" %in% stats) {
      tmp <- sapply( # mi x q x ntau
        1:ntau, \(itau) {
          asl(Phi_i,itau) %*% (
            diag(lambda[,itau],q,q) %*% asl(Phit_Phi,itau) +
              sigma2[itau] * diag(1,q,q))
        }, simplify = "array"
      )
      grad_Theta <- grad_Theta +  # p x ntau
        2 / n * sapply(
          1:ntau, \(itau) {
            as.matrix(Matrix::crossprod(Bi, asl(tmp,itau) %*% diag(lambda[,itau],q,q)))},
          simplify = "array") -
        2 / n * sapply(
          1:ntau, \(itau) {
            Bt_Yi %*% crossprod(Phit_Yi[,itau], diag(lambda[,itau],q,q))},
          simplify = "array")
      
      tr_PhiLamPhi <- sapply(1:ntau, \(itau) {
        sum(asl(Phi_i,itau)^2 %*% lambda[,itau])})  # ntau
      grad_zeta <- grad_zeta - sapply(1:ntau, \(itau) { # ntau
        sigma2[itau] / n * (
          Yit_Yi - mi * sigma2[itau] - tr_PhiLamPhi[itau])})
      for (k in 1:q) {
        grad_eta[k,] <- grad_eta[k,] - sapply(1:ntau, \(itau) {
          lambda[k,itau] / n * (
            Phit_Yi[k,itau]^2 - sigma2[itau] * sum(Phi_i[,k,itau]^2) -
              sum(Phit_Phi[,k,itau]^2 * lambda[,itau]))})
      }
    }
    
    if ("mse" %in% stats) {
      Lam_Phit_Phi <- sapply(1:ntau, \(itau) {  # q x q x ntau
        diag(lambda[,itau],q,q) %*% asl(Phit_Phi,itau)},
        simplify = "array")
      mse <- mse + 0.5 / n * sapply(1:ntau, \(itau) {
        Yit_Yi^2 + mi * sigma2[itau]^2 +
          sum(t(asl(Lam_Phit_Phi,itau)) * asl(Lam_Phit_Phi,itau)) -
          2 * sigma2[itau] * Yit_Yi -
          2 * sum(Phit_Yi[,itau]^2 * lambda[,itau]) +
          2 * sigma2[itau] * sum(asl(Phi_i,itau)^2 %*% lambda[,itau])
      })
    }
    
    if ("mse_offdiag" %in% stats) {
      Lam_Phit_Phi <- sapply(1:ntau, \(itau) {  # q x q x ntau
        diag(lambda[,itau],q,q) %*% asl(Phit_Phi,itau)},
        simplify = "array")
      diag_sse <- sapply(1:ntau, \(itau) {  # ntau
        Yit_Yi - sum(asl(Phi_i,itau)^2 %*% lambda[,itau])})
      mse_offdiag <- mse_offdiag +
        sapply(1:ntau, \(itau) {
          Yit_Yi^2 + sum(t(asl(Lam_Phit_Phi,itau)) * asl(Lam_Phit_Phi,itau)) -
            2 * sum(Phit_Yi[,itau]^2 * lambda[,itau]) - diag_sse[itau]})
    }
  }
  
  if ("grad" %in% stats) {
    if (!multiSmoothness) {
      grad_Theta <- grad_Theta[,,1]
      grad_eta <- grad_eta[,1]
    }
    out[['grad_Theta']] <- grad_Theta
    out[['grad_zeta']] <- grad_zeta
    out[['grad_eta']] <- grad_eta
  }
  if ("mse" %in% stats) {
    out[['mse']] <- mse
  }
  if ("mse_offdiag" %in% stats) {
    out[['mse_offdiag']] <- mse_offdiag
  }
  return(out)
}


pspc_smoothpen <- function(Theta, Omega) {
  q <- ncol(Theta)
  return(sum(sapply(1:q, \(k) sum(Theta[,k] * Omega %*% Theta[,k]))))
}

pspc_smoothpen_grad <- function(Theta, Omega) {
  return(2 * as.matrix(Omega %*% Theta))
}



