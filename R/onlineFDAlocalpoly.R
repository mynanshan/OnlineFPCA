#### Online local polynomials  #####

fpca.lpoly.online <- function(data_generator, n = 20, meanfun = TRUE,
                              evalArgs, streamArgs,
                              L2 = 3, L1 = L2, G = 0.5, R = 1,
                              optimal_C = TRUE, C0 = 1,
                              Mcl = 1, trueParams = NULL,
                              verbose = FALSE, period = 1) {
  ## online covariance function estimation
  EV1 <- evalArgs$EV1
  EV2 <- evalArgs$EV2
  t0 <- evalArgs$t0
  t1 <- evalArgs$t1
  eval_mu <- evalArgs$eval_mu
  eval_gam_vec <- evalArgs$eval_gam_vec
  eval_gam_mat <- evalArgs$eval_gam_mat
  Kmax <- streamArgs$Kmax
  # sds <- streamArgs$sds

  # # FIXME: select C
  # if (is.character(C)) {
  #   selected_C <- TRUE
  #   load(C) # 'sparse_cov_band_batch.Rdata'
  # } else {
  #   stopifnot(is.numeric(C) && C > 0)
  #   selected_C <- FALSE
  # }

  # initialize
  N <- 0
  mfull <- 0
  N_gam <- 0
  sigma_gam <- 0
  C_old <- 10
  sigma_gam3 <- 0
  sigma_gam5 <- 0

  # for mean est
  res_theta_mu <- list()
  res_theta_mu$centroids <- rep(0, L1)
  res_theta_mu$P <- array(0, dim = c(4, 4, EV1, L1))
  res_theta_mu$q <- array(0, dim = c(4, EV1, L1))
  res_sigma_mu1 <- list()
  res_sigma_mu1$centroids <- rep(0, L1)
  res_sigma_mu1$P <- array(0, dim = c(2, 2, EV1, L1))
  res_sigma_mu1$q <- array(0, dim = c(2, EV1, L1))
  res_sigma_mu2 <- res_sigma_mu1
  res_mu <- res_sigma_mu1
  # for cov est
  res_gam <- list()
  res_gam$centroids <- rep(0, L2)
  res_gam$P <- array(0, dim = c(3, 3, EV2^2, L2))
  res_gam$q <- array(0, dim = c(3, EV2^2, L2))
  # for C select
  res_sigma_gam1 <- list()
  res_sigma_gam1$centroids <- rep(0, 3)
  res_sigma_gam1$P <- array(0, dim = c(3, 3, EV2^2, 3))
  res_sigma_gam1$q <- array(0, dim = c(3, EV2^2, 3))
  res_sigma_gam2 <- res_sigma_gam1
  res_sigma_gam3 <- res_sigma_gam1
  res_theta_gam <- list()
  res_theta_gam$centroids <- rep(0, 3)
  res_theta_gam$P <- array(0, dim = c(10, 10, EV2^2, 3))
  res_theta_gam$q <- array(0, dim = c(10, EV2^2, 3))

  time <- c()
  # rss1 <- c()
  h1 <- c()
  theta <- c()
  sigma <- c()
  mus <- c()
  # rss2 <- c()
  vecs <- c()
  vals <- c()
  h2 <- c()
  Cs <- c()
  # time_band <- c()
  # C_old <- 1
  # dens<-list(); deris<-list()
  # s3<-c(); s5<-c()

  theta_mu <- 0
  sigma_mu <- 0
  h_old <- 1

  total_count <- 0

  for (K in 1:Kmax) {
    if (verbose) print(paste("K =", K))

    # generate data
    # set.seed(sds[K])
    dat <- data_generator(n, total_count)
    total_count <- total_count + n
    mk <- n
    njk1 <- dat$Lmi
    x <- unlist(dat$Lt)
    y <- unlist(dat$Ly)
    NK <- length(y)
    N <- N + NK
    mfull <- mfull + mk

    ### mean -------------------

    ### theta_mu
    h_theta_mu <- G * N^(-1 / 7)
    res_theta_mu <- online_LCub( # update sufficient statistics
      x, y, eval_mu, h_theta_mu, L1,
      res_theta_mu, N, NK, 1,
      is.1stBlock = (K == 1)
    )
    mu_sec_deri <- sapply(1:EV1, function(i) {
      2 * (solve(res_theta_mu$P[, , i, 1] + diag(1e-12, 4)) %*% matrix(res_theta_mu$q[, i, 1], 4, 1))[3]
    })
    theta_mu <- # optimal bandwidth
      mean(res_theta_mu$P[1, 1, 6:95, 1] * mu_sec_deri[6:95]^2) # 3117

    ### sigma_mu
    h_sigma_mu <- G * N^(-1 / 5)
    res_sigma_mu1 <- online_LL(x, y, eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1, is.1stBlock = (K == 1))
    mu <- sapply(1:EV1, function(i) {
      (solve(res_sigma_mu1$P[, , i, 1] + diag(1e-12, 2)) %*% matrix(res_sigma_mu1$q[, i, 1], 2, 1))[1]
    })
    mu_est <- approx(eval_mu, mu, xout = x, method = "linear")$y
    r <- (y - mu_est)^2
    res_sigma_mu2 <- online_LL(x, r, eval_mu, h_sigma_mu, L1, res_sigma_mu2, N, NK, 1, is.1stBlock = (K == 1))
    r <- sapply(1:EV1, function(i) {
      (solve(res_sigma_mu2$P[, , i, 1] + diag(1e-12, 2)) %*% matrix(res_sigma_mu2$q[, i, 1], 2, 1))[1]
    })
    sigma_mu <- (N - NK) / N * sigma_mu + NK / N * mean(r) # 0.87

    ### update mean
    h_mu <- min((15 * sigma_mu / theta_mu)^(1 / 5) * N^(-1 / 5), h_old)
    res_mu <- online_LL(x, y, eval_mu, h_mu, L1, res_mu, N, NK, 1, is.1stBlock = (K == 1))
    mu <- sapply(1:EV1, function(i) {
      (solve(res_mu$P[, , i, 1] + diag(1e-12, 2)) %*% matrix(res_mu$q[, i, 1], 2, 1))[1]
    })
    h_old <- h_mu

    # record
    # rss1 <- c(rss1, mean((mu - mu_true)^2))
    h1 <- c(h1, h_mu)
    theta <- c(theta, theta_mu)
    sigma <- c(sigma, sigma_mu)
    mus <- cbind(mus, mu)

    time.start <- Sys.time()

    ### covariance --------------------------
    mu <- mus[, K]
    mu_est <- approx(eval_mu, mu, xout = x, method = "linear")$y
    res_gam_data <- gene_gam_data(NK, njk1, mk, x, y, mu_est)
    u <- res_gam_data$u
    v <- res_gam_data$v
    NK_gam <- sum(njk1 * (njk1 - 1))
    N_gam <- N_gam + NK_gam

    ### theta_gam
    if (optimal_C) {
      h_theta_gam <- G * N_gam^(-1 / 8)
      res_theta_gam <- online_LCub(
        u, v, eval_gam_mat, h_theta_gam, 3,
        res_theta_gam, N_gam, NK_gam, 2
      )
      gam_sec_deri <- sapply(1:EV2^2, function(i) {
        2 * (sum((solve(res_theta_gam$P[, , i, 1] +
          diag(1e-12, 10)) %*%
          matrix(res_theta_gam$q[, i, 1], 10, 1))[c(4, 6)]))
      })
      gam_sec_deri <- matrix(gam_sec_deri, EV2, EV2)
      den <- matrix(res_theta_gam$P[1, 1, , 1], EV2, EV2)
      deri <- gam_sec_deri
      theta_gam <- mean(den[2:49, 2:49] * gam_sec_deri[2:49, 2:49]^2) # 242
      # dens<-c(dens,list(den))
      # deris<-c(deris,list(gam_sec_deri))

      ### sigma_gam
      h_sigma_gam <- G * N_gam^(-1 / 6)
      res_sigma_gam1 <- online_LL(
        u, v, eval_gam_mat, h_sigma_gam,
        3, res_sigma_gam1, N_gam, NK_gam, 2
      )
      gam <- sapply(1:EV2^2, function(i) {
        (solve(res_sigma_gam1$P[, , i, 1] + diag(1e-12, 3)) %*%
          matrix(res_sigma_gam1$q[, i, 1], 3, 1))[1]
      })
      r <- (v - sapply(1:NK_gam, function(i) {
        gam[which.min(abs(u[i, 1] - eval_gam_mat[, 1]) +
          abs(u[i, 2] - eval_gam_mat[, 2]))]
      }))^2
      r[r > (mean(r) + 5 * sqrt(var(r)))] <- mean(r)
      r[r < (mean(r) - 5 * sqrt(var(r)))] <- mean(r)
      res_sigma_gam2 <- online_LL(
        u, r, eval_gam_mat, h_sigma_gam,
        3, res_sigma_gam2, N_gam, NK_gam, 2
      )
      r <- sapply(1:EV2^2, function(i) {
        (solve(res_sigma_gam2$P[, , i, 1] + diag(1e-12, 3)) %*%
          matrix(res_sigma_gam2$q[, i, 1], 3, 1))[1]
      })
      r <- matrix(r, EV2, EV2)
      gam <- matrix(gam, EV2, EV2)
      sigma_eps <- sigma[K] - mean(diag(gam))
      V1 <- mean(r) + 2 * mean(diag(gam)) * sigma_eps + sigma_eps^2 - mean(gam^2)
      V1 <- V1 + mean(diag(r)) + 2 * mean(diag(gam)) * sigma_eps +
        sigma_eps^2 - mean(diag(gam)^2)
      sigma_gam5 <- (N_gam - NK_gam) / N_gam * sigma_gam5 +
        NK_gam / N_gam * V1 # 1.66
      # s5<-c(s5,sigma_gam5)

      ### C select
      s5 <- sigma_gam5
      th <- mean(den[3:47, 3:47] * deri[3:47, 3:47]^2)
      C <- (2 / 0.04 * 0.6^2 * s5 / th)^(1 / 6)
    } else {
      # disable optimal bandwidth selection to save some time
      C <- C0
    }

    Cs <- c(Cs, C)

    ### update covariance
    C_gam <- (1 - NK_gam / N_gam) * C_old + (NK_gam / N_gam) * C
    C_old <- C_gam
    h_gam <- C_gam * N_gam^(-1 / 6)
    res_gam <- online_LL(u, v, eval_gam_mat, h_gam, L2, res_gam, N_gam, NK_gam, 2)
    gam <- sapply(1:EV2^2, function(i) {
      (solve(res_gam$P[, , i, 1] + diag(1e-12, 3)) %*% matrix(res_gam$q[, i, 1], 3, 1))[1]
    })

    # principal component analysis
    vectors <- eigen(matrix(gam, EV2, EV2))$vectors
    values <- eigen(matrix(gam, EV2, EV2))$values
    mod <- sapply(1:EV2, function(i) sqrt(mean(vectors[, i]^2) * (t1 - t0)))
    vectors <- sapply(1:EV2, function(i) sign(values[i]) * vectors[, i] / mod[i])
    values <- values * mod^2
    N_pos_pc <- sum(values > 0)

    time.end <- Sys.time()

    # record
    if (K %% period == 0) {
      vecs <- cbind(vecs, as.vector(vectors))
      vals <- cbind(vals, values)
    }
    # rss2 <- c(rss2, mean((gam - gam_true)^2))
    h2 <- c(h2, h_gam)
    time <- c(time, difftime(time.end, time.start, units = "secs"))
  }

  result <- list(
    time = time, theta = theta, mus = mus, Cs = Cs,
    vecs = vecs, vals = vals, h1 = h1, h2 = h2
  )
  return(result)
}

#### Batch local polynomials  #####
# fpca.lpoly.batch <- function(Ly, Lt, meanfun = TRUE,
#                              evalArgs, streamArgs,
#                              L2 = 3, L1 = L2, G = 0.5, R = 1,
#                              C = 1,
#                              Mcl = 1, trueParams = NULL,
#                              verbose = FALSE) {
#   # not implemented as not needed here
# }


#### Computational functions ####

online_LCub <- function(x, y, eval, h, L, res_list, N, n, d,
                        is.1stBlock = FALSE) {
  ## sufficient statistics for the local linear smoother

  eta <- sapply(1:L, function(l) { # candidate bandwidths
    ((L - l + 1) / L)^(1 / (6 + d)) * h
  })

  if (is.1stBlock) { # connect bandwidth paths
    idx <- sapply(1:L, function(l) {
      which.min(abs(eta[l] - res_list$centroids))
    })
  } else {
    idx <- 1:L
  }

  # weighted average bandwidth
  res_list$centroids <- (res_list$centroids[idx] * (N - n) + eta * n) / N

  if (d == 1) {
    EV <- length(eval)

    for (l in 1:L) {
      Pnew <- array(0, dim = c(4, 4, EV))
      qnew <- matrix(0, 4, EV)

      for (i in 1:EV) {
        side <- cbind(1, x - eval[i], (x - eval[i])^2, (x - eval[i])^3)
        K_vec <- Epan((x - eval[i]) / eta[l]) / eta[l]
        for (nr in 1:4) {
          for (nc in 1:4) {
            Pnew[nr, nc, i] <- sum(K_vec * side[, nr] * side[, nc]) / n
          }
          qnew[nr, i] <- sum(K_vec * side[, nr] * y) / n
        }
      }

      res_list$P[, , , l] <- (res_list$P[, , , idx[l]] * (N - n) + Pnew * n) / N
      res_list$q[, , l] <- (res_list$q[, , idx[l]] * (N - n) + qnew * n) / N
    }
  } else {
    EV <- nrow(eval)

    for (l in 1:L) {
      Pnew <- array(0, dim = c(10, 10, EV))
      qnew <- matrix(0, 10, EV)
      for (i in 1:EV) {
        side <- cbind(
          1, x[, 1] - eval[i, 1], x[, 2] - eval[i, 2],
          (x[, 1] - eval[i, 1])^2,
          (x[, 1] - eval[i, 1]) * (x[, 2] - eval[i, 2]),
          (x[, 2] - eval[i, 2])^2,
          (x[, 1] - eval[i, 1])^3,
          (x[, 1] - eval[i, 1])^2 * (x[, 2] - eval[i, 2]),
          (x[, 1] - eval[i, 1]) * (x[, 2] - eval[i, 2])^2,
          (x[, 2] - eval[i, 2])^3
        ) # dim: N*10
        K_vec <- (Epan((x[, 1] - eval[i, 1]) / eta[l]) / eta[l]
          * Epan((x[, 2] - eval[i, 2]) / eta[l]) / eta[l])
        for (nr in 1:10) {
          for (nc in 1:10) {
            Pnew[nr, nc, i] <- sum(K_vec * side[, nr] * side[, nc]) / n
          }
          qnew[nr, i] <- sum(K_vec * side[, nr] * y) / n
        }
      }

      res_list$P[, , , l] <- (res_list$P[, , , idx[l]] * (N - n) + Pnew * n) / N
      res_list$q[, , l] <- (res_list$q[, , idx[l]] * (N - n) + qnew * n) / N
    }
  }

  return(res_list)
}


online_LL <- function(x, y, eval, h, L, res_list, N, n, d, is.1stBlock = FALSE) {
  eta <- sapply(1:L, function(l) {
    ((L - l + 1) / L)^(1 / (d + 4)) * h
  })

  if (is.1stBlock) {
    idx <- sapply(1:L, function(l) {
      which.min(abs(eta[l] - res_list$centroids))
    })
  } else {
    idx <- 1:L
  }

  res_list$centroids <- (res_list$centroids[idx] * (N - n) + eta * n) / N

  if (d == 1) {
    EV <- length(eval)

    for (l in 1:L) {
      Pnew <- array(0, dim = c(2, 2, EV))
      qnew <- matrix(0, 2, EV)

      for (i in 1:EV) {
        side <- cbind(1, x - eval[i])
        K_vec <- Epan((x - eval[i]) / eta[l]) / eta[l]
        Pnew[, , i] <- matrix(c(
          sum(K_vec * side[, 1]^2), sum(K_vec * side[, 1] * side[, 2]),
          sum(K_vec * side[, 2] * side[, 1]), sum(K_vec * side[, 2]^2)
        ), 2, 2) / n
        qnew[, i] <- matrix(c(
          sum(K_vec * side[, 1] * y), sum(K_vec * side[, 2] * y)
        ), 2, 1) / n
      }

      res_list$P[, , , l] <- (res_list$P[, , , idx[l]] * (N - n) + Pnew * n) / N
      res_list$q[, , l] <- (res_list$q[, , idx[l]] * (N - n) + qnew * n) / N
    }
  } else {
    EV <- nrow(eval)

    for (l in 1:L) {
      Pnew <- array(0, dim = c(3, 3, EV))
      qnew <- matrix(0, 3, EV)

      for (i in 1:EV) {
        side <- cbind(1, x[, 1] - eval[i, 1], x[, 2] - eval[i, 2])
        K_vec <- (Epan((x[, 1] - eval[i, 1]) / eta[l]) / eta[l]
          * Epan((x[, 2] - eval[i, 2]) / eta[l]) / eta[l])
        Pnew[, , i] <- matrix(c(
          sum(K_vec * side[, 1]^2), sum(K_vec * side[, 1] * side[, 2]), sum(K_vec * side[, 1] * side[, 3]),
          sum(K_vec * side[, 2] * side[, 1]), sum(K_vec * side[, 2]^2), sum(K_vec * side[, 2] * side[, 3]),
          sum(K_vec * side[, 3] * side[, 1]), sum(K_vec * side[, 3] * side[, 2]), sum(K_vec * side[, 3]^2)
        ), 3, 3) / n
        qnew[, i] <- matrix(c(
          sum(K_vec * side[, 1] * y), sum(K_vec * side[, 2] * y),
          sum(K_vec * side[, 3] * y)
        ), 3, 1) / n
      }

      res_list$P[, , , l] <- (res_list$P[, , , idx[l]] * (N - n) + Pnew * n) / N
      res_list$q[, , l] <- (res_list$q[, , idx[l]] * (N - n) + qnew * n) / N
    }
  }

  return(res_list)
}


batch_LCub <- function(x, y, eval, h, N, d) {
  if (d == 1) {
    EV <- length(eval)

    P <- array(0, dim = c(4, 4, EV))
    q <- matrix(0, 4, EV)

    for (i in 1:EV) {
      side <- cbind(1, x - eval[i], (x - eval[i])^2, (x - eval[i])^3)
      K_vec <- Epan((x - eval[i]) / h) / h
      for (nr in 1:4) {
        for (nc in 1:4) {
          P[nr, nc, i] <- sum(K_vec * side[, nr] * side[, nc]) / N
        }
        q[nr, i] <- sum(K_vec * side[, nr] * y) / N
      }
    }

    sec_deri <- sapply(1:EV, function(i) {
      2 * (solve(P[, , i] + diag(1e-12, 4)) %*% matrix(q[, i], 4, 1))[3]
    })

    index1 <- 6 * (EV == EV1) + 3 * (EV == EV2)
    index2 <- 95 * (EV == EV1) + 48 * (EV == EV2)

    theta <- mean(P[1, 1, index1:index2]
    * sec_deri[index1:index2]^2)
  } else {
    EV <- nrow(eval)

    P <- array(0, dim = c(10, 10, EV))
    q <- matrix(0, 10, EV)
    for (i in 1:EV) {
      side <- cbind(
        1, x[, 1] - eval[i, 1], x[, 2] - eval[i, 2],
        (x[, 1] - eval[i, 1])^2,
        (x[, 1] - eval[i, 1]) * (x[, 2] - eval[i, 2]),
        (x[, 2] - eval[i, 2])^2,
        (x[, 1] - eval[i, 1])^3,
        (x[, 1] - eval[i, 1])^2 * (x[, 2] - eval[i, 2]),
        (x[, 1] - eval[i, 1]) * (x[, 2] - eval[i, 2])^2,
        (x[, 2] - eval[i, 2])^3
      )
      K_vec <- (Epan((x[, 1] - eval[i, 1]) / h) / h
        * Epan((x[, 2] - eval[i, 2]) / h) / h)
      for (nr in 1:10) {
        for (nc in 1:10) {
          P[nr, nc, i] <- sum(K_vec * side[, nr] * side[, nc]) / N
        }
        q[nr, i] <- sum(K_vec * side[, nr] * y) / N
      }
    }

    sec_deri <- sapply(1:EV, function(i) {
      2 * (sum((solve(P[, , i] + diag(1e-12, 10)) %*% matrix(q[, i], 10, 1))[c(4, 6)]))
    })
    sec_deri <- matrix(sec_deri, EV2, EV2)
    den <- matrix(P[1, 1, ], EV2, EV2)
    theta <- mean(den[2:48, 2:48] * sec_deri[2:48, 2:48]^2) # 242

    return(theta)
  }
}

batch_LL <- function(x, y, eval, h, N, d) {
  if (d == 1) {
    EV <- length(eval)

    P <- array(0, dim = c(2, 2, EV))
    q <- matrix(0, 2, EV)

    for (i in 1:EV) {
      side <- cbind(1, x - eval[i])
      K_vec <- Epan((x - eval[i]) / h) / h
      P[, , i] <- matrix(c(
        sum(K_vec * side[, 1]^2), sum(K_vec * side[, 1] * side[, 2]),
        sum(K_vec * side[, 2] * side[, 1]), sum(K_vec * side[, 2]^2)
      ), 2, 2) / N
      q[, i] <- matrix(c(
        sum(K_vec * side[, 1] * y), sum(K_vec * side[, 2] * y)
      ), 2, 1) / N
    }

    est <- sapply(1:EV, function(i) {
      (solve(P[, , i] + diag(1e-12, 2)) %*% matrix(q[, i], 2, 1))[1]
    })
  } else {
    EV <- nrow(eval)

    P <- array(0, dim = c(3, 3, EV))
    q <- matrix(0, 3, EV)

    for (i in 1:EV) {
      side <- cbind(1, x[, 1] - eval[i, 1], x[, 2] - eval[i, 2])
      K_vec <- (Epan((x[, 1] - eval[i, 1]) / h) / h
        * Epan((x[, 2] - eval[i, 2]) / h) / h)
      P[, , i] <- matrix(c(
        sum(K_vec * side[, 1]^2), sum(K_vec * side[, 1] * side[, 2]), sum(K_vec * side[, 1] * side[, 3]),
        sum(K_vec * side[, 2] * side[, 1]), sum(K_vec * side[, 2]^2), sum(K_vec * side[, 2] * side[, 3]),
        sum(K_vec * side[, 3] * side[, 1]), sum(K_vec * side[, 3] * side[, 2]), sum(K_vec * side[, 3]^2)
      ), 3, 3) / N
      q[, i] <- matrix(c(
        sum(K_vec * side[, 1] * y), sum(K_vec * side[, 2] * y),
        sum(K_vec * side[, 3] * y)
      ), 3, 1) / N
    }

    est <- sapply(1:EV, function(i) {
      (solve(P[, , i] + diag(1e-12, 3)) %*% matrix(q[, i], 3, 1))[1]
    })
  }

  return(est)
}


# reform observations to two-dim coordinates and raw covariances
gene_gam_data <- function(NK, njk, mk, x, y, mu_est) {
  coor1 <- unlist(sapply(1:NK, function(i) {
    rep((1:NK)[i], njk[min(which(i <= cumsum(njk)))] - 1)
  }))
  II <- c(0, cumsum(njk))
  coor2 <- unlist(sapply(1:mk, function(i) {
    a <- (1:NK)[(II[i] + 1):II[i + 1]]
    b <- as.vector(sapply(1:length(a), function(j) a[-j]))
  }))
  rm(II)
  u <- t(sapply(1:length(coor1), function(i) c(x[coor1[i]], x[coor2[i]])))
  v <- sapply(1:length(coor1), function(i) {
    (y[coor1[i]] - mu_est[coor1[i]]) * (y[coor2[i]] - mu_est[coor2[i]])
  })
  rm(coor1, coor2)

  res <- list(u, v)
  names(res) <- c("u", "v")
  return(res)
}


# Epanechnikov Kernel
Epan <- function(z) {
  return(3 / 4 * (1 - z^2) * (abs(z) < 1))
}
