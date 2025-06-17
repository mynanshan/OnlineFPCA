# Regression-type FPCA for sparse and irregular data

library(fda)
library(Matrix)


fpca.reg <- function(
    Ly, Ltid, inits, meanfun = FALSE, npc=3,
    tau=1e-6, nu = 2, delta = 1,
    maxIter = 20, verbose = TRUE, sampleIter = 20,
    record_iterations = FALSE, record_period = 5,
    use_validation_set = FALSE, n_valid_ratio = 0.3,
    refine_alpha = FALSE) {
  
  # inits: Theta, alphamat
  
  tau_list = sort(tau)
  ntau = length(tau_list)
  
  N <- length(Ly)
  if (use_validation_set) {
    Ntest = round(N * n_valid_ratio)
    testIdx = sample(seq_len(N), Ntest)
    testLy = Ly[testIdx]
    testLtid = Ltid[testIdx]
    testLmi = sapply(testLy, length)
    Ly = Ly[-testIdx]
    Ltid = Ltid[-testIdx]
  }

  N <- length(Ly)
  Lmi <- sapply(Ly, length)
  m <- sum(Lmi)
  
  # arguments check
  stopifnot(!is.null(inits$Theta))
  stopifnot(is.numeric(maxIter))
  
  Theta <- inits$Theta  # p x q
  theta_mu <- inits$theta_mu
  if (!is.matrix(Theta) && is.vector(Theta)) Theta <- matrix(Theta, ncol=1)
  
  p <- nrow(Theta)
  q <- ncol(Theta)
  while (q < npc) {
    theta_append = runif(p, -1, 1)
    theta_append = manifold.Stiefel.retractSinglePC(theta_append, Theta, G)
    if (sum(theta_append^2) < .Machine$double.eps) next()
    Theta = cbind(Theta, theta_append)
    q = ncol(Theta)
  }
  # stopifnot(npc <= q)
  # q <- npc  # TODO: when q<npc, random initialize rest PCs
  stopifnot(q == ncol(Theta))
  stopifnot(N == length(Ltid))
  stopifnot(all(Lmi == sapply(Ltid, length)))
  if (is.null(theta_mu) || !meanfun) theta_mu <- rep(0, p)
  
  # useful invariants in global env
  # basis matrices: B
  # \int b(t) b(t)^\top dt: G
  # penalty matrix \int b''(t) b''(t)^\top dt: Omega
  
  nrecord = ceiling(maxIter / record_period)
  maxIter = nrecord * record_period
  if (record_iterations) objvals <- array(NA, dim = c(q, nrecord, ntau))
  if (record_iterations) Theta_history <- array(NA, dim = c(p, q, nrecord, ntau))
  Theta_est = array(NA, dim = c(p, q, ntau))
  alphamat <- matrix(NA, nrow=N, ncol=q)

  Lz = mapply(function(Yi, Ti) {
    Yi - as.matrix(B[Ti, , drop = F] %*% theta_mu)
  }, Ly, Ltid)
  if (use_validation_set) {
    testLz = mapply(function(Yi, Ti) {
      Yi - as.matrix(B[Ti, , drop = F] %*% theta_mu)
    }, testLy, testLtid)
    
    test_scores = matrix(NA, nrow = q, ncol = ntau)
  }

  conv_flag = matrix(0, nrow=q, ncol=ntau)
  nIter = matrix(0, nrow = q, ncol = ntau)

  gcv_scores = gcv_scores_alpha = gcv_scores_theta =
    matrix(NA, nrow = q, ncol = ntau)
  ## choice 1: based on fitting difficulty
  theta_gcv_factor1 = p / (N + p)  # larger dim of theta, higher difficulty for theta
  theta_gcv_factor2 = 1 - prod(Lmi^{-1/N})  # less N, higher mi, higher difficulty for theta
  theta_gcv_factor = sqrt(theta_gcv_factor1 * theta_gcv_factor2)
  ## choice 2: based on #{numbers it affects}
  # theta_gcv_factor = m / (N + m)
  ## choice 3: just root square them
  # theta_gcv_factor = 0.3
  ## choice 4: just use theta
  # theta_gcv_factor = 0.99
  tau_select = rep(NA, q)
  mse_all = matrix(NA, nrow=q, ncol=ntau)
  trSa_all = matrix(NA, nrow=q, ncol=ntau)
  trSth_all = matrix(NA, nrow=q, ncol=ntau)

  # main iteratios -----------------------------

  for (k in 1:q) {
    if (verbose) message("Update FPC ", k)

    thetak_init = Theta[, k]

    for (itau in seq_along(tau_list)) {
      
      tau = tau_list[itau]
      
      thetak = thetak_init
      phik_rough = sum(thetak * Omega %*% thetak)
      Lphik = lapply(Ltid, function(Ti) {
        as.vector(B[Ti, , drop = F] %*% thetak)
      })

      for (ii in 1:maxIter) {

        if (verbose && ii %% record_period == 0) message("Iter ", ii)

        thetak_old = thetak
        ir = (ii - 1) %/% record_period + 1

        # compute FPC scores
        alphak = fpca.reg.getAlphak(
          Lz, Lphik, phik_rough, Lmi, tau, nu, delta)

        # M: sum_i alpha_ik^2 B_i^T B_i / m_i
        # u: sum_i alpha_ik B_i^T z_{i,k-1} / m_i
        M <- matrix(0, nrow = p, ncol = p)
        u <- matrix(0, nrow = p, ncol = 1)
        for (i in 1:N) {
          mi <- Lmi[i]
          Bi <- B[Ltid[[i]], , drop = F] # mi x p
          Zi <- Lz[[i]]

          M <- M + alphak[i]^2 * as.matrix(crossprod(Bi)) / mi / N # p x p
          u <- u + alphak[i] * as.matrix(crossprod(Bi, Zi)) / mi / N # p x 1
        }

        # update fpc[k]
        thetak <- solve(
          as.matrix(M + tau * mean(alphak^2 * Lmi^(-1/(1+nu))) / N * Omega), u)

        # re-orthonormalize
        if (k > 1) {
          inprod_phik_phi <- as.matrix(crossprod(Theta[, 1:(k - 1), drop = F], G %*% thetak))
          thetak <- thetak - Theta[, 1:(k - 1), drop = F] %*% inprod_phik_phi
        }
        thetak <- thetak / sqrt(sum(thetak * G %*% thetak))

        # update fpc evaluations
        phik_rough = sum(thetak * Omega %*% thetak)
        Lphik = lapply(Ltid, function(Ti) {
          as.vector(B[Ti, , drop = F] %*% thetak)
        })

        # record
        Theta_est[, k, itau] = thetak
        if (record_iterations && (ii %% record_period == 0)) {
          Theta_history[, k, ir, itau] <- thetak
          objvals[k, ir, itau] = mapply(function(Zi, alpha_ik, phi_ik) {
            mean((Zi - alpha_ik * phi_ik)^2) # SSE / mi
          }, Lz, alphak, Lphik, SIMPLIFY = TRUE) |> mean() # SSE / mi / N
        }

        # check convergence
        if (max(abs(thetak - thetak_old)) / max(abs(thetak_old)) < 1e-4 ||
          max(abs(thetak - thetak_old)) < 1e-4) {
          break
        }
      } # end of main iterations

      if (verbose && ii == maxIter) message("FPC ", k, " has not converged")
      conv_flag[k, itau] = !(ii == maxIter) * 1
      nIter[k, itau] = ii
      
      # compute GCV related stats
      wsse = mapply(function(Zi, alpha_ik, phi_ik) {
        mean((Zi - alpha_ik * phi_ik)^2)  # SSE_i / mi
      }, Lz, alphak, Lphik, SIMPLIFY = TRUE) |> sum()
      num  = wsse / m
      # GCV: alphak as the variable
      phi_ik_norm = sapply(1:N, \(i) { sum(Lphik[[i]]^2)} )
      gamk = tau * phik_rough  # reformulated penalty coef
      tr_Salpha_div_m = sum(phi_ik_norm / (phi_ik_norm + gamk)) / m
      den_alpha = (1 - tr_Salpha_div_m)^2
      gcv_scores_alpha[k, itau] = num / den_alpha
      # GCV: thetak as the variable
      etak = tau * mean(Lmi^(-1 / (1 + nu)) * alphak^2)
      tr_Stheta_div_m = sum(diag(solve(M + etak * Omega, M))) / m
      den_theta = (1 - tr_Stheta_div_m)^2
      gcv_scores_theta[k, itau] = num / den_theta
      # aggregate two scores
      gcv_scores[k, itau] =
        gcv_scores_theta[k, itau]^theta_gcv_factor *
        gcv_scores_alpha[k, itau]^(1 - theta_gcv_factor)
      mse_all[k, itau] = wsse
      trSa_all[k, itau] = tr_Salpha_div_m
      trSth_all[k, itau] = tr_Stheta_div_m
      
      # compute test score
      if (use_validation_set) {
        testLphik = lapply(testLtid, function(Ti) {
          as.vector(B[Ti, , drop = F] %*% thetak)
        })
        alphak_test = fpca.reg.getAlphak(
          testLz, testLphik, phik_rough, testLmi, tau, nu, delta)
        test_scores[k, itau] = mapply(function(Zi, alpha_ik, phi_ik) {
          mean((Zi - alpha_ik * phi_ik)^2) # SSE / mi
        }, testLz, alphak_test, testLphik, SIMPLIFY = TRUE) |> mean()
      }

    } # end of all tau
    
    # parameter selection
    if (!use_validation_set) {
      itau_min = which.min(gcv_scores[k, ])
    } else {
      itau_min = which.min(test_scores[k, ])
    }
    tau_select[k] = tau_list[itau_min]
    Theta[, k] = Theta_est[, k, itau_min]

    # update stored stats
    phik_rough = sum(Theta[, k] * Omega %*% Theta[, k])
    Lphik = lapply(Ltid, function(Ti) {
      as.vector(B[Ti, , drop = F] %*% Theta[, k])
    })
    
    # finalize the FPC scores
    if (!refine_alpha) {
      alphamat[, k] = fpca.reg.getAlphak(
        Lz, Lphik, phik_rough, Lmi, tau_select[k], nu, delta)
    } else {
      alphamat[, k] = mapply(function(Zi, phi_ik) {
        fastmatrix::ridge(
          y ~ x - 1, data = data.frame(y = Lz[[i]], x = Lphik[[i]]))$coefficients
      }, Lz, Lphik)
    }
    # update Lz for the next fpc estimation
    Lz = mapply(function(Zi, alpha_ik, phi_ik) {
      Zi - alpha_ik * phi_ik
    }, Lz, alphamat[, k], Lphik)
    
    if (use_validation_set) {
      # update stats for the validation set
      testLphik = lapply(testLtid, function(Ti) {
        as.vector(B[Ti, , drop = F] %*% Theta[, k])
      })
      alpha_test = mapply(
        function(Zi, phi_ik, mi) {
          (mean(phi_ik^2) + tau * mi^(-1 / (1 + nu)) * phik_rough)^(-1) *
            mean(phi_ik * Zi)
        }, testLz, testLphik, testLmi
      )
      testLz = mapply(function(Zi, alpha_ik, phi_ik) {
        Zi - alpha_ik * phi_ik
      }, testLz, alpha_test, testLphik)
    }

  } # end of all FPC

  return(list(
    Theta = Theta, theta_mu = theta_mu, alphamat = alphamat,
    records = if (record_iterations) {
      list(Theta_history = Theta_history, objvals = objvals)
    } else {NULL},
    tau_list = tau_list, tau_select = tau_select,
    nu=nu, delta=delta,
    gcv_scores = gcv_scores, gcv_scores_alpha = gcv_scores_alpha,
    gcv_scores_theta = gcv_scores_theta,
    conv_flag = conv_flag, theta_gcv_factor = theta_gcv_factor,
    test_scores = if (use_validation_set) test_scores else {NULL},
    test_stats = list(
      mse_all = mse_all, trSa_all = trSa_all, trSth_all = trSth_all
    )
  ))
}


fpca.reg.getAlphak = function(
    Lz, Lphik, phik_rough, Lmi, tau, nu, delta) {
  alphak = mapply(
    function(Zi, phi_ik, mi) {
      # alpha_jk =
      # (phi_ik^T phi_ik + \tau m_i^{nu/(1+nu)} beta_k^T Omega beta_k)^{-1} *
      # phi_ik^T z_{i,k-1}
      (mean(phi_ik^2) +
         tau * mi^(-1 / (1 + nu)) * (phik_rough + delta))^(-1) *
        mean(phi_ik * Zi)
    }, Lz, Lphik, Lmi
  )
  return(alphak)
}


fpca.reg.fitNewdata = function(
    Ly, Ltid, Theta, theta_mu=NULL,
    tau=1e-6, nu=2, delta=2) {
  
  N = length(Ly)
  q = ncol(Theta)
  Lmi = sapply(Ly, length)
  
  if (!is.null(theta_mu)) {
    Lz = mapply(function(Yi, Ti) {
      Yi - as.matrix(B[Ti, , drop = F] %*% theta_mu)
    }, Ly, Ltid)
  } else {Lz = Ly}
  
  alphamat = matrix(nrow=N, ncol=q)
  mse = numeric(q)
  objval = numeric(q)

  for (k in 1:q) {
    thetak = Theta[,k]
    phik_rough = sum(thetak * Omega %*% thetak)
    Lphik = lapply(Ltid, function(Ti) {
      as.vector(B[Ti, , drop = F] %*% thetak)
    })
    
    alphamat[, k] = fpca.reg.getAlphak(
      Lz, Lphik, phik_rough, Lmi, tau[k], nu, delta)
    
    Lz = mapply(function(Zi, alpha_ik, phi_ik) {
      Zi - alpha_ik * phi_ik
    }, Lz, alphamat[, k], Lphik)
    
    mse[k] = mean(sapply(Lz, \(x) mean(x^2)))
    objval[k] = mse[k] +
      tau[k] * mean(alphamat[,k]^2 * Lmi^(-1 / (1 + nu))) * phik_rough
  }
  
  return(list(alphamat = alphamat, mse = mse, objval = objval))
}



