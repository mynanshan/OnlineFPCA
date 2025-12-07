library(Matrix)

fpca.sgd <- function(
  data_generator,
  obsGrid,
  inits,
  meanfun = FALSE,
  npc = 3,
  tau = NULL,
  tau.control = list(),
  nbatch = 1,
  stepsize = 1e-1,
  maxIter = 1000,
  nIter.constStepSize = 0,
  stepsize.decayrate = 0.5,
  stepsize.min = 0,
  period.decay = 1,
  nIter.slowerdecay = floor(0.5 * maxIter),
  stepsize.decayrate.slow = 0,
  dynlr = FALSE,
  dynlrCtrl = NULL,
  sgdtype = c(
    "sgd",
    "sgdm",
    "adagrad",
    "adam",
    "adam0",
    "adagrad2",
    "adam2",
    "rasa"
  ),
  beta1 = 0.9,
  beta2 = 0.99,
  adamw = TRUE,
  adam.rescale = FALSE,
  ada.start = 1,
  adareset = 200,
  adareset.end = Inf,
  asgd.start = 1,
  asgd.reset = 200,
  asgd.reset.end = Inf,
  asgd.end = Inf,
  weight = c("obs", "subj"),
  fpcCI = FALSE,
  nIter.1stTune = floor(0.1 * maxIter),
  nIter.lastTune = maxIter,
  period.tune = 100,
  vcrit.tune = c("ewmabv", "av"),
  ewmabv.beta = NULL,
  abv_aofv_w = 0.5,
  nIter.tauNoIncrease = floor(0.6 * maxIter),
  period.message = 100,
  recordParams = TRUE,
  period.record = 20,
  period.time = period.record,
  verbose = TRUE,
  test.mode = FALSE
) {
  cat(">>>>>>>>>>> Start OnlineFPCA <<<<<<<<<<<\n")
  # Set params and initializations -----------------------------------------

  inits <- setParams.inits(inits, meanfun)
  optns.tau <- setParams.tau(tau, tau.control)
  tau <- sort(optns.tau$tau, decreasing = TRUE)
  taurange0 <- ifelse(length(tau) > 1, diff(range(log(tau))), 1)
  tau.control <- optns.tau$tau.control
  ntau <- tau.control$ntau
  adatau <- tau.control$adatau
  nselect <- tau.control$nselect
  nchild <- allocate_nchild(ntau, nselect)
  cat("Dynamic tuning:", adatau, "\n")
  if (adatau) {
    tau.history <- tau
    tau.selectId <- c()
  } else {
    tau.history <- NULL
    tau.selectId <- NULL
  }
  ntune <- floor((maxIter - nIter.1stTune) / period.tune) + 1
  if (is.null(ewmabv.beta)) {
    if (ntune < 5) {
      ewmabv.beta <- 0.5
    } else {
      ewmabv.beta.candidates <- seq(0.1, 0.95, 0.05)
      valid.ewmabv.beta.id <- which(
        ewmabv.beta.candidates^(ntune - 1) < 0.5 / ntune^1.2
      )
      if (length(valid.ewmabv.beta.id) == 0) {
        idx <- 1
        warning(
          "Improper choice of ewmabv.beta. The tuning period is possibly too large."
        )
      }
      idx <- max(valid.ewmabv.beta.id)
      if (is.null(ewmabv.beta) || ewmabv.beta > ewmabv.beta.candidates[idx]) {
        ewmabv.beta <- ewmabv.beta.candidates[idx]
      }
    }
  } else {
    stopifnot(ewmabv.beta <= 1 && ewmabv.beta >= 0)
  }
  vcrit.tune <- match.arg(vcrit.tune)
  vcrit.history <- list(
    bv = array(dim = c(ntau, ntune)),
    av = array(dim = c(ntau, ntune)),
    ewmabv = array(dim = c(ntau, ntune))
  )
  vcrit.history$bv.avg <- array(dim = c(ntau, ntune))
  vcrit.history$av.avg <- array(dim = c(ntau, ntune))
  vcrit.history$ewmabv.avg <- array(dim = c(ntau, ntune))
  time.history <- array(0, dim = c(ntau, maxIter))

  iter.select <- 1
  stopifnot(stepsize > 0)
  stopifnot(stepsize.decayrate > 0 && stepsize.decayrate <= 1)
  stopifnot(stepsize.decayrate.slow >= 0 && stepsize.decayrate.slow <= 1)
  if (nIter.slowerdecay < nIter.constStepSize) {
    nIter.slowerdecay <- nIter.constStepSize
  }

  p <- nrow(inits$Theta)
  q <- ncol(inits$Theta)
  while (q < npc) {
    theta_append = runif(p, -1, 1)
    theta_append = manifold.Stiefel.normalize(
      theta_append,
      inits$Theta,
      G
    )
    if (sum(theta_append^2) < .Machine$double.eps) {
      next()
    }
    inits$Theta = cbind(inits$Theta, theta_append)
    q = ncol(inits$Theta)
  }
  Theta <- array(inits$Theta, dim = c(p, q, ntau))
  lambda <- array(inits$lambda, dim = c(q, ntau))
  sigma2 <- array(inits$sigma2, dim = c(ntau))
  theta_mu <- inits$theta_mu

  # Iteratively averaged estimates
  Theta.avg <- array(NA, dim = c(p, q, ntau))
  lambda.avg <- array(NA, dim = c(q, ntau))
  sigma2.avg <- array(NA, dim = c(ntau))
  asgd_counter <- 0

  if (recordParams) {
    nrecord <- floor(maxIter / period.record) + 1
    params.history <- list(
      Theta = array(dim = c(p, q, ntau, nrecord)),
      lambda = array(dim = c(q, ntau, nrecord)),
      sigma2 = array(dim = c(ntau, nrecord))
    )
    params.history$Theta[,,, 1] <- Theta
    params.history$lambda[,, 1] <- lambda
    params.history$sigma2[, 1] <- sigma2
    params.history$Theta.avg <- params.history$Theta
    params.history$lambda.avg <- params.history$lambda
    params.history$sigma2.avg <- params.history$sigma2
  } else {
    params.history <- NULL
  }
  iter.params <- 1

  # TODO: implement online mean estimation
  if (!is.null(theta_mu)) {
    stop("NOT IMPLEMENTED")
  }

  # storing Euclidean Hessians for vectorized parameters
  cat("Online CI:", fpcCI, "\n")
  vh <- if (!fpcCI) {
    NULL
  } else {
    list(
      H_Theta_Theta = array(0, dim=c(p*q, p*q)),
      H_eta_Theta = array(0, dim=c(q, p*q)),
      H_zeta_Theta = array(0, dim=c(p*q)),
      H_Theta_eta = matrix(0, p*q, q),
      H_eta_eta = matrix(0, q, q),
      H_zeta_eta = numeric(q),
      H_Theta_zeta = numeric(p*q),
      H_eta_zeta = numeric(q),
      H_zeta_zeta = numeric(1)
    )
  }
  # storing Euclidean gradients
  ge <- if (!fpcCI) {
    NULL
  } else {
    list(
      grad_Theta = matrix(0, p, q),
      grad_eta = numeric(q),
      grad_zeta = numeric(1)
    )
  }
  phi.var <- if (!fpcCI) {
    NULL
  } else {
    array(0, dim = c(nrow(B), q))
  }
  CIs <- if (!fpcCI) {
    NULL
  } else {
    array(0, dim = c(nrow(B), q, 3, nrecord))
  }
  if (fpcCI) {
    class(vh) <- "vecHess"
    attr(CIs, "iters") <- c()
  }

  # SGD and Adaptive SGD settings
  sgdtype <- match.arg(sgdtype)
  ada_stats <- init_ada_stats(p, q, ntau, sgdtype, adamw)
  init_grad <- FALSE
  if (all(c("grad_Theta", "grad_eta", "grad_zeta") %in% names(inits))) {
    init_grad <- TRUE
    for (itau in seq_len(itau)) {
      ada_stats <- update_ada_stats(
        ada_stats,
        inits$grad_Theta,
        inits$grad_eta,
        inits$grad_zeta,
        asl(Theta, l),
        NULL,
        itau = itau,
        i = 0,
        sgdtype = sgdtype,
        adamw = adamw
      )
    }
  }
  ada_counter <- as.integer(0 + init_grad)
  # exponential moving average of updating directions' norms
  ewadn <- numeric(ntau)

  # step size settings
  if (dynlr) {
    if (is.null(dynlrCtrl)) dynlrCtrl <- list()
    if (is.null(dynlrCtrl$niter)) {
      dynlrCtrl$niter <- 500
    } else {
      stopifnot(dynlrCtrl$niter > 0)
    }
    if (is.null(dynlrCtrl$reset)) {
      dynlrCtrl$reset <- 100
    } else {
      stopifnot(dynlrCtrl$reset > 0)
    }
    if (is.null(dynlrCtrl$refdn)) {
      dynlrCtrl$refdn <- 0.2
    } else {
      stopifnot(dynlrCtrl$refdn > 0)
    }
    if (is.null(dynlrCtrl$w)) {
      dynlrCtrl$w <- 0.99
    } else {
      stopifnot(dynlrCtrl$w > 0 && dynlrCtrl$w < 1)
    }
    direc_Th_norms <- rep(1, q)
    dynlr_counter <- 0
  }

  # ABV accumulators for dynamic tuning
  av <- numeric(ntau) # averaged one-function validation score
  ewmabv <- numeric(ntau) # averaged block validation score
  curr_bv <- numeric(ntau) # current block validation score
  av.avg <- numeric(ntau)
  ewmabv.avg <- numeric(ntau)
  curr_bv.avg <- numeric(ntau)

  # other settings
  weight <- match.arg(weight)

  # Main Iterations --------------------------------------------------------

  # record number of samples processed
  total_count <- 0

  # record some debugging info
  check <- list(
    grad_norms = c(),
    direc_norms = c(),
    stepsize = c(),
    dTheta = c(),
    mhat_norms = c(),
    vhat_sqrt = c(),
    m_norms = c(),
    v_sqrt = c()
  )

  decay_counter <- 0
  for (i in 1:maxIter) {
    if (verbose && (i %% period.message == 0)) {
      message("- Iteration ", i)
    }

    if (i == ada.start && !stringr::str_detect(sgdtype, "sgd")) {
      cat("Start adagrad-type optimization.\n")
      if (dynlr) {
        stepsize <- dynlrCtrl$refdn
        cat("Changing the step size to", stepsize, "\n")
      }
    }

    # Stepsize decaying factor
    if (i == nIter.constStepSize + 1) {
      decay_counter <- 0
    } else if (i == nIter.slowerdecay + 1) {
      ii <- max(decay_counter - 1, 1) %/% period.decay + 1
      # lr_new / ii^alpha = lr / ii^alpha.slow
      # ==> lr_new = lr / ii^alpha.slow * ii^alpha
      stepsize <- stepsize * ii^(stepsize.decayrate - stepsize.decayrate.slow) 
      # decay_counter <- 0
    }
    decay_counter <- decay_counter + 1
    decay <- 1 / (max(decay_counter - 1, 1) %/% period.decay + 1)^
      if (i <= nIter.constStepSize) {
        0
      } else if (i <= nIter.slowerdecay) {
        stepsize.decayrate.slow
      } else {
        stepsize.decayrate
      }
    curr_stepsize <- stepsize * decay
    if (curr_stepsize < stepsize.min) {
      curr_stepsize <- stepsize.min
    }

    if (asgd_counter == 0 && i >= asgd.start) {
      # if asgd_counter == 0, initialize avg estimates
      Theta.avg <- Theta
      lambda.avg <- lambda
      sigma2.avg <- sigma2
    }

    # receive new data
    obs <- data_generator(nbatch, total_count)
    Ltid <- obs$Ltid
    Ly <- obs$Ly
    total_count <- total_count + nbatch

    # calculate validation score
    for (l in 1:ntau) {
      # if (i > nIter.lastTune && l > 1) break
      if (i < asgd.start) {
        lik <- objfun(
          Ly,
          Ltid,
          asl(Theta, l),
          lambda[, l],
          sigma2[l],
          theta_mu,
          tau[l],
          stats = "loss"
        )$lik
        av[l] <- av[l] * (i - 1) / i + lik * 1 / i
        curr_bv[l] <- curr_bv[l] + lik
      } else {
        lik.avg <- objfun(
          Ly,
          Ltid,
          asl(Theta.avg, l),
          lambda.avg[, l],
          sigma2.avg[l],
          theta_mu,
          tau[l],
          stats = "loss"
        )$lik
        av.avg[l] <- av.avg[l] * (i - 1) / i + lik.avg * 1 / i
        curr_bv.avg[l] <- curr_bv.avg[l] + lik.avg
      }
    }

    if (
      i > 1 && i >= nIter.1stTune && (i - nIter.1stTune) %% period.tune == 0
    ) {
      itune <- floor((i - nIter.1stTune) / period.tune) + 1
      if (itune == 1) {
        curr_bv <- curr_bv / (i * nbatch)
      } else {
        curr_bv <- curr_bv / (period.tune * nbatch)
      }
      if (i < asgd.start) {
        vcrit.history$bv[, itune] <- curr_bv
        ewmabv <- ewmabv * ewmabv.beta + curr_bv * (1 - ewmabv.beta)
        curr_bv <- numeric(ntau)
        vcrit.history$av[, itune] <- av
        vcrit.history$ewmabv[, itune] <- ewmabv
      } else {
        curr_bv.avg <- curr_bv.avg / (period.tune * nbatch)
        vcrit.history$bv.avg[, itune] <- curr_bv.avg
        ewmabv.avg <- ewmabv.avg * ewmabv.beta + curr_bv.avg * (1 - ewmabv.beta)
        curr_bv.avg <- numeric(ntau)
        vcrit.history$av.avg[, itune] <- av.avg
        vcrit.history$ewmabv.avg[, itune] <- ewmabv.avg
      }

      # online selection of tau
      if (adatau && i <= nIter.lastTune) {
        iter.select <- c(iter.select, i)
        if (i < asgd.start) {
          vscore <- eval(parse(text = vcrit.tune))
        } else {
          vscore <- eval(parse(text = paste0(vcrit.tune, ".avg")))
        }
        lstar <- order(vscore)[1:nselect]
        parentId <- rep(lstar, nchild)
        tau.selectId <- unname(cbind(tau.selectId, parentId))
        ada_stats <- select_ada_stats(ada_stats, parentId)
        av <- av[parentId]
        ewmabv <- ewmabv[parentId]
        Theta <- Theta[,, parentId, drop = F]
        lambda <- lambda[, parentId, drop = F]
        sigma2 <- sigma2[parentId]
        tau.control <- setParams.tau(
          tau,
          tau.control,
          delta.min = taurange0 / max(1, ntau - 2) / (itune^0.7)
        )$tau.control
        tau <- explore_newtau(
          tau[lstar],
          nchild,
          tau.control$delta,
          i > nIter.tauNoIncrease
        )
        tau.history <- unname(cbind(tau.history, tau))
        if (i >= asgd.start) {
          Theta.avg <- Theta.avg[,, parentId, drop = F]
          lambda.avg <- lambda.avg[, parentId, drop = F]
          sigma2.avg <- sigma2.avg[parentId]
          av.avg <- av.avg[parentId]
          ewmabv.avg <- ewmabv.avg[parentId]
        }
      }
    }

    for (l in 1:ntau) {
      time.start <- Sys.time()

      # compute gradients
      objective <- objfun(
        Ly,
        Ltid,
        asl(Theta, l),
        lambda[, l],
        sigma2[l],
        theta_mu,
        tau[l],
        stats = "grad",
        weight = weight
      )
      grad_Theta <- objective$grad_Theta
      grad_Theta <- as.matrix(manifold.Stiefel.project(
        grad_Theta,
        asl(Theta, l),
        G
      ))
      grad_eta <- objective$grad_eta
      grad_zeta <- objective$grad_zeta

      # update adagrad accumulators
      ada_stats <- update_ada_stats(
        ada_stats,
        grad_Theta,
        grad_eta,
        grad_zeta,
        asl(Theta, l),
        NULL,
        itau = l,
        i = ada_counter,
        sgdtype = sgdtype,
        adamw = adamw
      )

      # formulate updating directions
      direc <- get_ada_direc(
        grad_Theta,
        grad_eta,
        grad_zeta,
        ada_stats,
        l,
        sgdtype = ifelse(i >= ada.start, sgdtype, "sgd"),
        adamw = adamw
      )
      if (sgdtype != "sgd" && i >= ada.start) {
        # if not sgd, the direction may be not a tangent element
        direc[['Theta']] <- as.matrix(manifold.Stiefel.project(
          direc[['Theta']],
          asl(Theta, l),
          G
        ))
      }
      # TODO: testing ewadn. This looks like a dirty fix. Can we get rid of it?
      if (
        stringr::str_detect(sgdtype, "adam") && i >= ada.start  # Adam always need a direction standardization
      ) {
        mn_dTheta <- sum(direc[['Theta']] * G %*% direc[['Theta']]) / q
        ewadn[l] <- ifelse(
          ada_counter == 0 && i >= ada.start,
          mn_dTheta,
          ewadn[l] * 0.9 + 0.1 * mn_dTheta
        )
        direc[['Theta']] <- direc[['Theta']] / sqrt(ewadn[l])
      }

      # store information for CIs
      # CI's accumulation should be subject to the restart of ASGD
      if (fpcCI && i >= asgd.start && i <= asgd.end && l == 1) {
        # Only use l==1 for CI computation
        ii <- asgd_counter  # num of mini-batches used for CIs
        # accumulate Euclidean gradients
        for (gname in names(ge)) {
          ge[[gname]] <- ge[[gname]] * ii / (ii + 1) +
            objective[[gname]] / (ii + 1)
        }
        # accumulate Euclidean Hessian
        vhi <- computeVecHess(
          Ly, Ltid,
          asl(Theta, 1), lambda[, 1], sigma2[1],
          theta_mu, tau[1]
        )
        for (vname in names(vh)) {
          vh[[vname]] <- vh[[vname]] * ii / (ii + 1) +
            vhi[[vname]] / (ii + 1)
        }
        curr_Theta <- asl(Theta.avg, 1)
        Delta_Th <- manifold.Stiefel.project(
          ge$grad_Theta, curr_Theta, G
        )
        # TODO: compute inverse Hessian, switch to a CG solver
        vecGrad <- c(Delta_Th, ge$grad_eta, ge$grad_zeta)
        opres <- suppressWarnings(
          gslnls::gsl_nls(
            fn = \(x) vecHess_apply_wrapper(x, vh, curr_Theta, grad_Theta),
            y = vecGrad,
            start = vecGrad
          )
        )
        vsol <- vecProj(opres$m$getPars(), curr_Theta)
        V <- matrix(vsol[1:(p*q)], p, q)
        JV <- dV2phi(V, matrix(0,p,q), curr_Theta)
        curr_phi_var <- JV^2
        phi.var <- phi.var * ii / (ii + 1) +
          curr_phi_var / (ii + 1)
        # update CI estimates
        if (i %% period.record == 0) {
          attr(CIs, "iters") <- c(attr(CIs, "iters"), i)
          idx <- i %/% period.record + 1
          curr_sd <- sqrt(phi.var) / sqrt(ii + 1)
          curr_Phi <- as.matrix(B %*% curr_Theta)
          CIs[,,1,idx] <- curr_Phi - 1.96 * curr_sd
          CIs[,,2,idx] <- curr_Phi
          CIs[,,3,idx] <- curr_Phi + 1.96 * curr_sd
          if (verbose) cat("CI computed at iteration", i, "\n")
        }
      }

      # Parameter updates
      Theta_l <- asl(Theta, l)
      Theta[,, l] <- manifold.Stiefel.retract(
        -curr_stepsize * direc[['Theta']],
        Theta_l, G
      )
      lambda[, l] <- lambda[, l] * exp(-curr_stepsize * direc[['other']][1:q])
      sigma2[l] <- sigma2[l] * exp(-curr_stepsize * direc[['other']][q + 1])

      if (i >= asgd.start && i <= asgd.end) {
        invRetr_Theta <- manifold.Stiefel.invRetract(
          asl(Theta, l),
          asl(Theta.avg, l),
          G
        )
        Theta.avg[,, l] <- manifold.Stiefel.retract(
          invRetr_Theta / (asgd_counter + 1),
          asl(Theta.avg, l),
          G
        )
        lambda.avg[, l] <- lambda.avg[, l]^(asgd_counter / (asgd_counter + 1)) *
          lambda[, l]^(1 / (asgd_counter + 1))
        sigma2.avg[l] <- sigma2.avg[l]^(asgd_counter / (asgd_counter + 1)) *
          sigma2[l]^(1 / (asgd_counter + 1))
      }

      if (dynlr && l == 1 && i <= dynlrCtrl$niter) {
        dThL2 <- colSums(direc[['Theta']] * G %*% direc[['Theta']])
        direc_Th_norms <- if (dynlr_counter == 0) {
          dThL2
        } else {
          dThL2 * (1 - dynlrCtrl$w) + direc_Th_norms * dynlrCtrl$w
          # dThL2 / (dynlr_counter + 1) + direc_Th_norms * dynlr_counter / (dynlr_counter + 1)
        }
        dynlr_counter <- dynlr_counter + 1
        if (i %% dynlrCtrl$reset == 0) {
          stepsize <- dynlrCtrl$refdn / (decay * sqrt(mean(direc_Th_norms)))
          if (verbose) {
            cat("Step size updated to", stepsize, "\n")
          }
          dynlr_counter <- 0
          direc_Th_norms <- rep(0, q)
        }
      }

      # In-Fitting informations. May be changed to a callback 
      if (verbose && l == 1 && (i %% period.message == 0)) {
        direc_norms <- sqrt(c(
          Theta = colSums(direc[['Theta']] * G %*% direc[['Theta']]) |> sqrt(),
          eta = (direc[['other']][1:q]) |> abs(),
          zeta = direc[['other']][q + 1] |> abs()
        ))
        print(direc_norms * curr_stepsize)
      }

      if (l == 1) {
        direc_norms <- sqrt(c(
          Theta = colSums(direc[['Theta']] * G %*% direc[['Theta']]) |> sqrt(),
          eta = (direc[['other']][1:q]) |> abs(),
          zeta = direc[['other']][q + 1] |> abs()
        ))
        check$grad_norms <- rbind(check$grad_norms, sqrt(colSums(grad_Theta * G %*% grad_Theta)))
        check$direc_norms <- rbind(check$direc_norms, sqrt(colSums(direc[['Theta']] * G %*% direc[['Theta']])))
        check$stepsize <- c(check$stepsize, curr_stepsize)
        check$dTheta <- rbind(check$dTheta, sqrt(colMeans(as.matrix(B %*% (Theta[,, l] - Theta_l))^2)))
        # check$mhat_norms <- rbind(check$m_norms, sqrt(colSums(ada_stats$mhat$Theta[,,l] * G %*% ada_stats$mhat$Theta[,,l])))
        # check$vhat_sqrt <- rbind(check$vhat_sqrt, sqrt(diag(ada_stats$vhat$Theta[,,l])))
        # check$m_norms <- rbind(check$m_norms, sqrt(colSums(ada_stats$m$Theta[,,l] * G %*% ada_stats$m$Theta[,,l])))
        # check$v_sqrt <- rbind(check$v_sqrt, sqrt(diag(ada_stats$v$Theta[,,l])))
      }

      # record time
      time.end <- Sys.time()
      time.history[l, i] <- time.history[l, i] +
        difftime(time.end, time.start, units = 'secs')
    } # end l in 1:ntau

    # update adaptive gradient counter
    # if (i >= ada.start) {
      if (i <= adareset.end && i %% adareset == 0) {
        ada_counter <- 0
        ada_stats <- init_ada_stats(p, q, ntau, sgdtype, adamw)
      } else {
        ada_counter <- ada_counter + 1
      }
    # }

    # update averaging counter
    if (i >= asgd.start && i <= asgd.end) {
      if (i <= asgd.reset.end && i %% asgd.reset == 0) {
        asgd_counter <- 0
        # push back the averaged estimates
        Theta <- Theta.avg
        lambda <- lambda.avg
        sigma2 <- sigma2.avg
      } else {
        asgd_counter <- asgd_counter + 1
      }
    }

    # record parameters
    if (recordParams && (i %% period.record == 0)) {
      iter.params <- c(iter.params, i)
      idx <- i %/% period.record + 1
      params.history$Theta[,,, idx] <- Theta
      params.history$lambda[,, idx] <- lambda
      params.history$sigma2[, idx] <- sigma2
      params.history$Theta.avg[,,, idx] <- Theta.avg
      params.history$lambda.avg[,, idx] <- lambda.avg
      params.history$sigma2.avg[, idx] <- sigma2.avg
    }

  } # end of iteration i

  # Final selection
  minId <- if (i >= asgd.start && i <= asgd.end) {
    which.min(ewmabv.avg)
  } else {
    which.min(ewmabv)
  }
  tau.min <- tau[minId]

  out <- list(
    Theta = Theta,
    lambda = lambda,
    sigma2 = sigma2,
    Theta.avg = Theta.avg,
    lambda.avg = lambda.avg,
    sigma2.avg = sigma2.avg,
    theta_mu = theta_mu,
    CI = CIs,
    tau = tau,
    tau.min = tau.min,
    av = av,
    ewmabv = ewmabv,
    av.avg = av.avg,
    ewmabv.avg = ewmabv.avg,
    tau.select = list(
      tau.history = tau.history,
      tau.selectId = tau.selectId,
      iter.select = iter.select
    ),
    params.history = list(params = params.history, iter.params = iter.params),
    vcrit.history = vcrit.history,
    time.history = time.history,
    stepsize = stepsize,
    stepsize.final = stepsize * decay,
    ewmabv.beta = ewmabv.beta,
    check = check
  )

  return(out)
}

objfun <- function(
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3,
  stats = c("all", "loss", "grad"),
  weight = c("obs", "subj")
) {
  # global constants: B, Omega
  stats <- match.arg(stats)
  weight <- match.arg(weight)
  p <- ncol(B)
  q <- length(lambda)
  n <- length(Ly)
  stopifnot(ncol(Theta) == q)
  stopifnot(nrow(Theta) == p)
  stopifnot(length(Ltid) == n)
  if (length(tau) == 1 && q > 1) {
    tau <- rep(tau, q)
  }
  stopifnot(length(tau) == q)
  if (stats == "all" || stats == "loss") {
    fval <- 0
  }
  if (stats == "all" || stats == "grad") {
    grad_Theta <- matrix(0, p, q)
    grad_eta <- numeric(q)
    grad_zeta <- 0
  }

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

    if (stats == "all" || stats == "loss") {
      logdetQ <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)
      fval <- fval +
        1 /
          n *
          (1 /
            sigma2 *
            (sum(Yi^2) - sum(PhiT_Yi * invQ %*% PhiT_Yi)) +
            (mi - q) * log(sigma2) +
            sum(log(lambda)) +
            logdetQ)
    }
    if (stats == "all" || stats == "grad") {
      Phit_Phi <- crossprod(Phi)
      invQ_PhiT_Yi <- invQ %*% PhiT_Yi
      grad_Theta <- grad_Theta +
        2 /
          n *
          as.matrix(
            crossprod(Bi, Phi) %*%
              (1 / sigma2 * invQ_PhiT_Yi %*% t(invQ_PhiT_Yi) + invQ) -
              1 / sigma2 * crossprod(Bi, Yi) %*% t(invQ_PhiT_Yi)
          )

      grad_eta <- grad_eta +
        1 /
          n *
          (-invQ_PhiT_Yi^2 / lambda + rep(1, q) - sigma2 * diag(invQ) / lambda)

      grad_zeta <- grad_zeta +
        1 /
          n *
          (-sum(Yi^2) /
            sigma2 +
            sum(invQ_PhiT_Yi * PhiT_Yi) / sigma2 +
            sum(invQ_PhiT_Yi^2 / lambda) +
            mi -
            q +
            sigma2 * sum(diag(invQ) / lambda))

      if (weight == "subj") {
        grad_Theta <- grad_Theta / mi
        grad_eta <- grad_eta / mi
        grad_zeta <- grad_zeta / mi
      }
    }
  }
  out <- list()
  if (stats == "all" || stats == "loss") {
    out[["fval"]] <- fval + sum((Theta * Omega %*% Theta) %*% tau)
    out[["lik"]] <- fval # unpenalized likelihood
  }
  if (stats == "all" || stats == "grad") {
    out[["grad_Theta"]] <- invG %*% (
      grad_Theta +
        2 * matrix(rep(tau, each = p), nrow = p) *
        as.matrix(Omega %*% Theta)
    ) |>
      as.matrix()
    out[["grad_eta"]] <- as.numeric(grad_eta)
    out[["grad_zeta"]] <- as.numeric(grad_zeta)
  }
  return(out)
}


# Support for Hessians and CIs -------------------------------------------

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
  Theta = NULL, grad_Theta = NULL,
  manifold = TRUE
) {
  stopifnot(is(vh, "vecHess"))
  if (manifold) {
    stopifnot(!is.null(Theta))
    stopifnot(!is.null(grad_Theta))
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

# TODO: Implement of the inverse Hessian
# Currently, we are using a dirty method: use the nonlinear ls
# also to solve the inverser.
# Ideally, we should have a conjugate gradient method

vecProj <- function(x, Theta) {
  Delta <- matrix(x[1:(p*q)], p, q)
  Delta <- manifold.Stiefel.project(Delta, Theta, G)
  c(Delta, tail(x, q + 1))
}

vecHess_apply_wrapper <- function(x, vh, Theta, grad_Theta) {
  x <- vecProj(x, Theta)
  Del_Th <- matrix(x[1:(p*q)], p ,q)
  del_eta <- x[(p*q+1):(p*q+q)]
  del_zeta <- x[(p*q+q+1):(p*q+q+1)]
  vha <- vecHess_apply(
    vh, Del_Th, del_eta, del_zeta,
    Theta, grad_Theta, TRUE
  )
  return(c(
    vha$H_Theta_Psi_Delta,
    vha$H_eta_Psi_Delta,
    vha$H_zeta_Psi_Delta
  ))
}

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

# Computational helpers --------------------------------------------------
# Some computational helpers
invG_mul <- function(X) {
  backsolve(GR, forwardsolve(t(GR), X))
}


# Raw mean and covj --------------------------------------------------------

get.sample.mean <- function(Ly, Ltid) {
  m <- length(unique(unlist(Ltid)))
  mean.hat <- numeric(m)
  mean.count <- numeric(m)
  for (i in seq_along(Ly)) {
    mean.hat[Ltid[[i]]] = mean.hat[Ltid[[i]]] + Ly[[i]]
    mean.count[Ltid[[i]]] = mean.count[Ltid[[i]]] + 1
  }
  mean.hat = mean.hat / mean.count
  mean.hat
}

get.sample.cov <- function(Ly, Ltid) {
  m <- length(unique(unlist(Ltid)))
  cov.hat <- matrix(0, nrow = m, ncol = m)
  cov.count <- matrix(0, nrow = m, ncol = m)
  for (i in seq_along(Ly)) {
    cov.hat[Ltid[[i]], Ltid[[i]]] = cov.hat[Ltid[[i]], Ltid[[i]]] +
      outer(Ly[[i]], Ly[[i]])
    cov.count[Ltid[[i]], Ltid[[i]]] = cov.count[Ltid[[i]], Ltid[[i]]] + 1
  }
  cov.hat = cov.hat / cov.count
  cov.hat
}


# AdaGrad-related helpers ------------------------------------------------

init_ada_stats <- function(
  p, q, ntau,
  sgdtype = c(
    "sgd",
    "sgdm",
    "adagrad",
    "adam",
    "adam0",
    "adagrad2",
    "adam2",
    "rasa"
  ),
  adamw = TRUE
) {
  sgdtype <- match.arg(sgdtype)
  ada_stats <- NULL
  if (sgdtype == "sgdm" || stringr::str_detect(sgdtype, "adam")) {
    ada_stats[['m']] <- list(
      Theta = array(0, dim = c(p, q, ntau)),
      other = matrix(0, nrow = q + 1, ncol = ntau)
    )
    if (adamw && stringr::str_detect(sgdtype, "adam")) {
      ada_stats[['mhat']] <- list(
        Theta = array(0, dim = c(p, q, ntau)),
        other = matrix(0, nrow = q + 1, ncol = ntau)
      )
    }
  }
  if (sgdtype %in% c("adagrad", "adam")) {
    ada_stats[['v']] <- list(
      Theta = array(0, dim = c(q, q, ntau)),
      other = array(0, dim = c(q + 1, q + 1, ntau))
    )
    if (adamw && sgdtype == "adam") {
      ada_stats[['vhat']] <- list(
        Theta = array(0, dim = c(q, q, ntau)),
        other = array(0, dim = c(q + 1, q + 1, ntau))
      )
    }
  } else if (sgdtype %in% c("adam0")) {
    ada_stats[['v']] <- list(
      Theta = array(0, dim = c(p, q, ntau)),
      other = array(0, dim = c(q + 1, ntau))
    )
    if (adamw) {
      ada_stats[['vhat']] <- list(
        Theta = array(0, dim = c(p, q, ntau)),
        other = array(0, dim = c(q + 1, ntau))
      )
    }
  } else if (sgdtype %in% c("adagrad2", "adam2")) {
    ada_stats[['v']] <- list(
      Theta.l = array(0, dim = c(p, ntau)),
      Theta.r = array(0, dim = c(q, q, ntau)),
      other = array(0, dim = c(q + 1, q + 1, ntau))
    )
    if (adamw && sgdtype == "adam2") {
      ada_stats[['vhat']] <- list(
        Theta.l = array(0, dim = c(p, ntau)),
        Theta.r = array(0, dim = c(q, q, ntau)),
        other = array(0, dim = c(q + 1, q + 1, ntau))
      )
    }
  } else if (sgdtype %in% c("rasa")) {
    ada_stats[['v']] <- list(
      Theta.l = array(0, dim = c(p, ntau)),
      Theta.r = array(0, dim = c(q, ntau)),
      other = array(0, dim = c(q + 1, ntau))
    )
    ada_stats[['vhat']] <- list(
      Theta.l = array(0, dim = c(p, ntau)),
      Theta.r = array(0, dim = c(q, ntau)),
      other = array(0, dim = c(q + 1, ntau))
    )
  }
  return(ada_stats)
}

update_ada_stats <- function(
  ada_stats,
  grad_Theta,
  grad_eta,
  grad_zeta,
  Theta,
  Theta_old = NULL,
  itau = 1,
  i = 0,
  sgdtype = c(
    "sgd",
    "sgdm",
    "adagrad",
    "adam",
    "adam0",
    "adagrad2",
    "adam2",
    "rasa"
  ),
  adamw = TRUE,
  beta1 = 0.9,
  beta2 = 0.99
) {
  sgdtype <- match.arg(sgdtype)
  p <- nrow(grad_Theta)
  q <- ncol(grad_Theta)
  stopifnot(length(grad_eta) == q)
  stopifnot(length(grad_zeta) == 1)

  # if (i == 0 && !stringr::str_detect(sgdtype, "adam")) {
  if (i == 0) {
    # If no gradient info is accummulated yet,
    # no averaging happens
    beta1 <- 0
    beta2 <- 0
  }

  # update momentum terms
  if (sgdtype == "sgdm" || stringr::str_detect(sgdtype, "adam")) {
    m_Theta <- ada_stats[['m']][['Theta']][,, itau]
    trs_m_Theta <- manifold.Stiefel.transport(m_Theta, Theta, NULL, G)
    m_Theta <- trs_m_Theta * beta1 + grad_Theta * (1 - beta1)
    m_other <- ada_stats[['m']][['other']][, itau]
    m_other <- m_other * beta1 + c(grad_eta, grad_zeta) * (1 - beta1)
    ada_stats[['m']][['Theta']][,, itau] <- m_Theta
    ada_stats[['m']][['other']][, itau] <- m_other
    if (stringr::str_detect(sgdtype, "adam") && adamw) {
      ada_stats[['mhat']][['Theta']][,, itau] <- m_Theta /
        ifelse(i > 0, 1 - beta1^i, 1)
      ada_stats[['mhat']][['other']][, itau] <- m_other /
        ifelse(i > 0, 1 - beta1^i, 1)
    }
  }

  # update variance terms for Theta
  if (sgdtype %in% c("adagrad", "adam")) {
    v_Theta <- ada_stats[['v']][['Theta']][,, itau]
    v_Theta_new <- as.matrix(t(grad_Theta) %*% G %*% grad_Theta)
    if (sgdtype == "adagrad") {
      ada_stats[['v']][['Theta']][,, itau] <-
        v_Theta * i / (i + 1) + v_Theta_new / (i + 1)
    } else {
      v_Theta <- v_Theta * beta2 + v_Theta_new * (1 - beta2)
      ada_stats[['v']][['Theta']][,, itau] <- v_Theta
      if (adamw) {
        ada_stats[['vhat']][['Theta']][,, itau] <- v_Theta /
          ifelse(i > 0, 1 - beta2^i, 1)
      }
    }
  } else if (sgdtype %in% c("adagrad2", "adam2")) {
    grad_Theta_til <- as.matrix(GR %*% grad_Theta)
    v_Theta_l <- ada_stats[['v']][['Theta.l']][, itau]
    v_Theta_r <- ada_stats[['v']][['Theta.r']][,, itau]
    v_Theta_l_new <- rowMeans(grad_Theta_til^2)     # p
    v_Theta_r_new <- crossprod(grad_Theta_til) / p  # q x q
    if (sgdtype == "adagrad2") {
      v_Theta_l <- v_Theta_l * i / (i + 1) + v_Theta_l_new / (i + 1)
      v_Theta_r <- v_Theta_r * i / (i + 1) + v_Theta_r_new / (i + 1)
    } else {
      v_Theta_l <- v_Theta_l * beta2 + v_Theta_l_new * (1 - beta2)
      v_Theta_r <- v_Theta_r * beta2 + v_Theta_r_new * (1 - beta2)
    }
    ada_stats[['v']][['Theta.l']][, itau] <- v_Theta_l
    ada_stats[['v']][['Theta.r']][,, itau] <- v_Theta_r
    if (adamw && sgdtype == "adam2") {
      ada_stats[['vhat']][['Theta.l']][, itau] <- v_Theta_l /
          ifelse(i > 0, 1 - beta2^i, 1)
      ada_stats[['vhat']][['Theta.r']][,, itau] <- v_Theta_r /
          ifelse(i > 0, 1 - beta2^i, 1)
    }
  } else if (sgdtype == "rasa") {
    grad_Theta_til <- as.matrix(GR %*% grad_Theta)
    v_Theta_l <-
      ada_stats[['v']][['Theta.l']][, itau] * beta2 +
      rowMeans(grad_Theta_til^2) * (1 - beta2)
    v_Theta_r <-
      ada_stats[['v']][['Theta.r']][, itau] * beta2 +
      colMeans(grad_Theta_til^2) * (1 - beta2)
    ada_stats[['v']][['Theta.l']][, itau] <- v_Theta_l
    ada_stats[['v']][['Theta.r']][, itau] <- v_Theta_r
    ada_stats[['vhat']][['Theta.l']][, itau] <-
      pmax(v_Theta_l, ada_stats[['vhat']][['Theta.l']][, itau])
    ada_stats[['vhat']][['Theta.r']][, itau] <-
      pmax(v_Theta_r, ada_stats[['vhat']][['Theta.r']][, itau])
  } else if (sgdtype == "adam0") {
    v_Theta <- ada_stats[['v']][['Theta']][,, itau]
    v_Theta_new <- grad_Theta^2
    v_Theta <- v_Theta * beta2 + v_Theta_new * (1 - beta2)
    ada_stats[['v']][['Theta']][,, itau] <- v_Theta
    if (adamw) {
      ada_stats[['vhat']][['Theta']][,, itau] <- v_Theta /
        ifelse(i > 0, 1 - beta2^i, 1)
    }
  }

  # update variance terms for eta and zeta
  if (sgdtype %in% c("adagrad", "adam", "adagrad2", "adam2")) {
    v_other <- ada_stats[['v']][['other']][,, itau]
    v_other_new <- outer(c(grad_eta, grad_zeta), c(grad_eta, grad_zeta))
    if (stringr::str_detect(sgdtype, "adagrad")) {
      ada_stats[['v']][['other']][,, itau] <-
        v_other * i / (i + 1) + v_other_new / (i + 1)
    } else {
      v_other <- v_other * beta2 + v_other_new * (1 - beta2)
      ada_stats[['v']][['other']][,, itau] <- v_other
      if (adamw) {
        ada_stats[['vhat']][['other']][,, itau] <- v_other /
          ifelse(i > 0, 1 - beta2^i, 1)
      }
    }
  } else if (sgdtype %in% c("rasa")) {
    v_other <- 
      ada_stats[['v']][['other']][, itau] * beta2 +
      c(grad_eta, grad_zeta)^2 * (1 - beta2)
    ada_stats[['v']][['other']][, itau] <- v_other
    ada_stats[['vhat']][['other']][, itau] <-
      pmax(v_other, ada_stats[['vhat']][['other']][, itau])
  } else if (sgdtype == "adam0") {
    v_other <- ada_stats[['v']][['other']][, itau]
    v_other_new <- c(grad_eta, grad_zeta)^2
    v_other <- v_other * beta2 + v_other_new * (1 - beta2)
    ada_stats[['v']][['other']][, itau] <- v_other
    if (adamw) {
      ada_stats[['vhat']][['other']][, itau] <- v_other /
        ifelse(i > 0, 1 - beta2^i, 1)
    }
  }

  ada_stats
}

safe_sqrtinv <- function(A, eps = 1e-6, eps2 = 1) {
  A <- 0.5 * (A + t(A))
  d <- nrow(A)
  eigvals <- eigen(A)$values
  if (any(eigvals < max(eigvals) * eps)) {
    eps <- eps2
  }
  pracma::sqrtm(A + diag(eps, d, d))$Binv
}

get_ada_direc <- function(
  grad_Theta,
  grad_eta,
  grad_zeta,
  ada_stats,
  itau,
  sgdtype = c(
    "sgd",
    "sgdm",
    "adagrad",
    "adam",
    "adam0",
    "adagrad2",
    "adam2",
    "rasa"
  ),
  adamw = TRUE
) {
  sgdtype = match.arg(sgdtype)
  direc <- list()
  if (sgdtype == "sgd") {
    direc[['Theta']] <- grad_Theta
    direc[['other']] <- c(grad_eta, grad_zeta)
  } else if (sgdtype == "sgdm") {
    direc[['Theta']] <- ada_stats[['m']][['Theta']][,, itau]
    direc[['other']] <- ada_stats[['m']][['other']][, itau]
  } else if (sgdtype == "adagrad") {
    # direc[['Theta']] <- grad_Theta %*%
    #   safe_sqrtinv(ada_stats[['v']][['Theta']][,, itau])
    direc[['Theta']] <- sweep(grad_Theta, 2,
      (diag(ada_stats[['v']][['Theta']][,, itau]) + 1e-8)^(-0.5), "*")
    # direc[['other']] <- 
    #   safe_sqrtinv(ada_stats[['v']][['other']][,, itau]) %*%
    #     c(grad_eta, grad_zeta)
    direc[['other']] <- c(grad_eta, grad_zeta) *
      (diag(ada_stats[['v']][['other']][,, itau]) + 1e-8)^(-0.5)
  } else if (sgdtype == "adam") {
    direc[['Theta']] <- if (adamw) {
      # ada_stats[['mhat']][['Theta']][,, itau] %*%
      #   safe_sqrtinv(ada_stats[['vhat']][['Theta']][,, itau])
      sweep(ada_stats[['mhat']][['Theta']][,, itau], 2,
        (diag(ada_stats[['vhat']][['Theta']][,, itau]) + 1e-8)^(-0.5), "*")
    } else {
      # ada_stats[['m']][['Theta']][,, itau] %*%
      #   safe_sqrtinv(ada_stats[['v']][['Theta']][,, itau])
      sweep(ada_stats[['m']][['Theta']][,, itau], 2,
        (diag(ada_stats[['v']][['Theta']][,, itau]) + 1e-8)^(-0.5), "*")
    }
    direc[['other']] <- if (adamw) {
      # safe_sqrtinv(ada_stats[['vhat']][['other']][,, itau]) %*%
      #   ada_stats[['mhat']][['other']][, itau]
      ada_stats[['mhat']][['other']][, itau] *
        (diag(ada_stats[['vhat']][['other']][,, itau]) + 1e-8)^(-0.5)
    } else {
      # safe_sqrtinv(ada_stats[['v']][['other']][,, itau]) %*%
      #   ada_stats[['m']][['other']][, itau]
      ada_stats[['m']][['other']][, itau] *
        (diag(ada_stats[['v']][['other']][,, itau]) + 1e-8)^(-0.5)
    }
  } else if (sgdtype == "adagrad2") {
    grad_Theta_til <- as.matrix(GR %*% grad_Theta)
    v_Theta_l <- ada_stats[['v']][['Theta.l']][, itau]
    v_Theta_r <- ada_stats[['v']][['Theta.r']][,, itau]
    v_other <- ada_stats[['v']][['other']][,, itau]
    # TODO: it is still not clear how to give 
    # norm-1 directions here. Think about it.
    # vsqrtinv_Theta_l <-
    #   (v_Theta_l * q / sum(diag(v_Theta_r)) + 1e-8)^(-0.5)
    # # vsqrtinv_Theta_r <-
    # #   safe_sqrtinv(v_Theta_r * p / sum(diag(v_Theta_r)))
    # vsqrtinv_Theta_r <- safe_sqrtinv(v_Theta_r)
    vsqrtinv_Theta_l <-
      (v_Theta_l * q / sum(diag(v_Theta_r)) + 1e-8)^(-0.5)
    vsqrtinv_Theta_r <- (diag(v_Theta_r) + 1e-8)^(-0.5)
    # direc[['Theta']] <- as.matrix(backsolve(
    #   GR,
    #   grad_Theta_til |> 
    #     sweep(1, vsqrtinv_Theta_l, "*") |> 
    #     (\(X) X %*% vsqrtinv_Theta_r)()
    # ))
    direc[['Theta']] <- as.matrix(backsolve(
      GR,
      grad_Theta_til |> 
        sweep(1, vsqrtinv_Theta_l, "*") |> 
        sweep(2, vsqrtinv_Theta_r, "*")
    ))
    # vsqrtinv_other <- safe_sqrtinv(v_other)
    # direc[['other']] <- vsqrtinv_other %*% c(grad_eta, grad_zeta)
    direc[['other']] <- c(grad_eta, grad_zeta) / sqrt(diag(v_other) + 1e-8)
  } else if (sgdtype == "adam2") {
    m_Theta <- if (!adamw) {
      ada_stats[['m']][['Theta']][,, itau]
    } else {
      ada_stats[['mhat']][['Theta']][,, itau]
    }
    m_Theta <- as.matrix(GR %*% m_Theta)
    m_other <- if (!adamw) {
      ada_stats[['m']][['other']][, itau]
    } else {
      ada_stats[['mhat']][['other']][, itau]
    }
    v_Theta_l <- ada_stats[['v']][['Theta.l']][, itau]
    v_Theta_r <- ada_stats[['v']][['Theta.r']][,, itau]
    v_other <- ada_stats[['v']][['other']][,, itau]
    vsqrtinv_Theta_l <-
      (v_Theta_l * q / sum(diag(v_Theta_r)) + 1e-8)^(-0.5)
    # vsqrtinv_Theta_r <-
    #   safe_sqrtinv(v_Theta_r * p / sum(diag(v_Theta_r)))
    vsqrtinv_Theta_r <- safe_sqrtinv(v_Theta_r)
    direc[['Theta']] <- as.matrix(backsolve(
      GR,
      m_Theta |> 
        sweep(1, vsqrtinv_Theta_l, "*") |> 
        (\(X) X %*% vsqrtinv_Theta_r)()
    ))
    vsqrtinv_other <- safe_sqrtinv(v_other)
    direc[['other']] <- vsqrtinv_other %*% c(grad_eta, grad_zeta)
  } else if (sgdtype == "rasa") {
    v_Theta_l <- ada_stats[['v']][['Theta.l']][, itau]
    v_Theta_r <- ada_stats[['v']][['Theta.r']][, itau]
    v_other <- ada_stats[['v']][['other']][, itau]
    v4rt_Theta_l <- (v_Theta_l + 1e-8)^(-0.25)
    v4rt_Theta_r <- (v_Theta_r + 1e-8)^(-0.25)
    grad_Theta_til <- as.matrix(GR %*% grad_Theta)
    direc[['Theta']] <- as.matrix(solve(
      GR,
      grad_Theta_til |> 
        sweep(1, v4rt_Theta_l, "*") |> 
        sweep(2, v4rt_Theta_r, "*")
    ))
    vsqrtinv_other <- (v_other + 1e-8)^(-0.5)
    direc[['other']] <- vsqrtinv_other * c(grad_eta, grad_zeta)
  } else if (sgdtype == "adam0") {
    direc[['Theta']] <- if (adamw) {
      ada_stats[['mhat']][['Theta']][,, itau] /
        sqrt(ada_stats[['vhat']][['Theta']][,, itau] + 1e-8)
    } else {
      ada_stats[['m']][['Theta']][,, itau] /
        sqrt(ada_stats[['v']][['Theta']][,, itau] + 1e-8)
    }
    direc[['other']] <- if (adamw) {
      ada_stats[['mhat']][['other']][, itau] /
        sqrt(ada_stats[['vhat']][['other']][, itau] + 1e-8)
    } else {
      ada_stats[['m']][['other']][, itau] /
        sqrt(ada_stats[['v']][['other']][, itau] + 1e-8)
    }
  }
  direc[['Theta']] <- as.matrix(direc[['Theta']])
  direc[['other']] <- as.vector(direc[['other']])
  direc
}

arrange_last_dim <- function(A, indices) {
  dims <- dim(A)
  ndim <- length(dims)
  last_dim <- dims[ndim]
  stopifnot(all(is.integer(indices)) && length(indices) == last_dim)
  array(matrix(A, ncol = last_dim)[, indices], dim = dims)
}

select_ada_stats <- function(ada_stats, parentId) {
  if (is.null(ada_stats)) {
    return(ada_stats)
  }
  lapply(
    ada_stats,
    \(a) {
      if (is.list(a)) {
        lapply(a, arrange_last_dim, indices = parentId)
      } else if (is.numeric(a)) {
        arrange_last_dim(a, indices = parentId)
      } else {
        stop("unknown type.")
      }
    }
  )
}

# Parameter validation ---------------------------------------------------

setParams.inits <- function(inits, meanfun = FALSE) {
  if (!("Theta" %in% names(inits))) {
    stop("Lack an initialization for Theta")
  }
  inits$Theta <- manifold.Stiefel.retract(inits$Theta, NULL, G)
  if (!("lambda" %in% names(inits))) {
    stop("Lack an initialization for lambda")
  }
  if (!("sigma2" %in% names(inits))) {
    stop("Lack an initialization for sigma2")
  }
  if (meanfun) {
    if (!("theta_mu" %in% names(inits))) {
      stop("Lack an initialization for theta_mu")
    }
  } else {
    inits$theta_mu <- NULL
  }
  return(inits)
}


# Dynamic tuning helpers -------------------------------------------------

setParams.tau <- function(
  tau = NULL, tau.control = list(),
  delta.min = NULL,
  k = 1,
  rate = 1
) {
  if (!is.null(tau)) {
    tau.control$ntau <- length(tau)
  } else {
    if (!"ntau" %in% names(tau.control)) {
      tau.control$ntau <- 9
    }
    if (!"maxtau" %in% names(tau.control)) {
      tau.control$maxtau <- 1e2
    }
    if (!"mintau" %in% names(tau.control)) {
      tau.control$mintau <- 1e-6
    }
    tau <- with(
      tau.control,
      exp(seq(log(mintau), log(maxtau), length.out = ntau))
    )
  }
  if (!"adatau" %in% names(tau.control)) {
    tau.control$adatau <- TRUE
  }
  if (tau.control$adatau) {
    if (!"nselect" %in% names(tau.control)) {
      tau.control$nselect <- 3
    }
    if (tau.control$nselect > tau.control$ntau) {
      tau.control$nselect <- tau.control$ntau
    }
    if (tau.control$nselect > tau.control$ntau / 2) {
      warning(
        "nselect > ntau/2 is not effective to explore the parameter space."
      )
    }
    # delta.min.default <- 0.3
    delta.min.default <- 0.1
    if (is.null(delta.min)) {
      delta.min <- delta.min.default
    } else {
      delta.min <- max(delta.min, delta.min.default)
    }
    if (length(tau) > 1) {
      # tau.control$delta <- max(delta.min, mean(diff(sort(log(tau)))) / 3)
      tau.control$delta <- max(delta.min, log(10) / k^rate)
    } else {
      tau.control$delta <- 0
    }
  }
  return(list(tau = tau, tau.control = tau.control))
}


allocate_nchild <- function(ntau, nselect, rule = "equal") {
  if (rule == "equal") {
    nchild <- rep(floor(ntau / nselect), nselect)
    k <- ntau - sum(nchild)
    if (k > 0) nchild[1:k] <- nchild[1:k] + 1
  }
  return(nchild)
}

explore_newtau <- function(tau, nchild, delta, noincrease = FALSE) {
  newtau <- c()
  for (l in seq_along(tau)) {
    if (noincrease) {
      newtau <- c(
        newtau,
        tau[l] * exp(-seq(0, 2 * delta, length.out = nchild[l]))
      )
    } else {
      newtau <- c(
        newtau,
        tau[l] * exp(-seq(-delta, delta, length.out = nchild[l]))
      )
    }
  }
  return(newtau)
}


plot.tau_path <- function(
  tau.history,
  tau.selectId,
  col.selected = "royalblue3",
  col.discarded = "gray50",
  lty.selected = 1,
  lty.discarded = 2,
  lwd.selected = 1.5,
  lwd.discarded = 1
) {
  d <- ncol(tau.history) - 1
  w <- nrow(tau.history)
  maxlogtau <- max(log10(tau.history))
  minlogtau <- min(log10(tau.history))
  uylim <- maxlogtau + 0.05 * (maxlogtau - minlogtau)
  lylim <- minlogtau - 0.05 * (maxlogtau - minlogtau)
  plot(
    x = NULL,
    xlim = c(0, d),
    ylim = c(lylim, uylim),
    xlab = expression(Number ~ of ~ Data ~ Blocks),
    ylab = expression(log[10](tau))
  )
  for (i in d:1) {
    taui <- tau.history[, i + 1]
    for (j in 1:length(taui)) {
      if (i < d && (!j %in% tau.selectId[, i + 1])) {
        curr_col <- col.discarded
        curr_lty <- lty.discarded
        curr_lwd <- lwd.discarded
      } else {
        curr_col <- col.selected
        curr_lty <- lty.selected
        curr_lwd <- lwd.selected
      }
      lines(
        x = c(i - 1, i),
        type = "b",
        col = curr_col,
        lty = curr_lty,
        lwd = curr_lwd,
        y = c(log10(tau.history[tau.selectId[j, i], i]), log10(taui[j]))
      )
    }
  }
  tau1 <- tau.history[, 1]
  tau1.selected <- tau1[tau.selectId[, 1]]
  tau1.discarded <- tau1[-tau.selectId[, 1]]
  points(
    x = rep(0, length(tau1.selected)),
    y = log10(tau1.selected),
    col = col.selected,
    lty = lty.selected,
    lwd = lwd.selected
  )
  points(
    x = rep(0, length(tau1.discarded)),
    y = log10(tau1.discarded),
    col = col.discarded,
    lty = lty.discarded,
    lwd = lwd.discarded
  )
}


plot.tau_path2 <- function(
  tau.history,
  tau.selectId,
  final.id = 1,
  nbatch_per_block = 1,
  col.selected = "royalblue3",
  col.discarded = "gray60",
  lty.selected = 1,
  lty.discarded = 2,
  lwd.selected = 1.5,
  lwd.discarded = 1,
  plot.title = NULL
) {
  d <- ncol(tau.history) - 1
  w <- nrow(tau.history)
  maxlogtau <- max(log10(tau.history))
  minlogtau <- min(log10(tau.history))
  uylim <- maxlogtau + 0.05 * (maxlogtau - minlogtau)
  lylim <- minlogtau - 0.05 * (maxlogtau - minlogtau)
  plot(
    x = NULL,
    xlim = c(0, d * nbatch_per_block),
    ylim = c(lylim, uylim),
    xlab = expression(Iteration ~ k),
    ylab = expression(log[10](tau[k])),
    main = plot.title
  )
  # draw all paths in grey
  for (i in d:1) {
    taui <- tau.history[, i + 1]
    for (j in 1:length(taui)) {
      lines(
        x = c(i - 1, i) * nbatch_per_block,
        type = "l",
        col = col.discarded,
        lty = lty.discarded,
        lwd = lwd.discarded,
        y = c(log10(tau.history[tau.selectId[j, i], i]), log10(taui[j]))
      )
    }
    points(
      rep(i, w) * nbatch_per_block,
      log10(taui),
      pch = 21,
      col = col.discarded,
      lwd = lwd.discarded,
      bg = "white"
    )
  }
  points(
    rep(0, w) * nbatch_per_block,
    log10(tau.history[, 1]),
    pch = 21,
    col = col.discarded,
    lwd = lwd.discarded,
    bg = "white"
  )
  # draw the selected path in blue
  tau_path <- extract_tau_path(tau.history, tau.selectId, final.id)
  for (i in d:1) {
    lines(
      x = c(i - 1, i) * nbatch_per_block,
      type = "l",
      col = col.selected,
      lty = lty.selected,
      lwd = lwd.selected,
      y = log10(c(tau_path$tau_path[i], tau_path$tau_path[i + 1]))
    )
    points(
      i * nbatch_per_block,
      log10(tau_path$tau_path[i + 1]),
      pch = 21,
      col = col.selected,
      lwd = lwd.selected,
      bg = "white"
    )
  }
  points(
    0 * nbatch_per_block,
    log10(tau_path$tau_path[1]),
    pch = 21,
    col = col.selected,
    lwd = lwd.selected,
    bg = "white"
  )
}


extract_tau_path <- function(tau.history, tau.selectId, idx = NULL) {
  d <- ncol(tau.history) - 1
  w <- nrow(tau.history)
  if (is.null(idx)) {
    idx <- 1
  }
  tau_path_id <- idx
  tau_path <- tau.history[idx, d + 1]
  for (i in d:1) {
    idx <- tau.selectId[idx, i]
    tau_path_id <- c(idx, tau_path_id)
    tau_path <- c(tau.history[idx, i], tau_path)
  }
  return(list(tau_path = tau_path, tau_path_id = tau_path_id))
}


# RMSE helper ------------------------------------------------------------

rmse_phi <- function(ThetaEst, ThetaTrue, B) {
  # ThetaEst: (p, q, L)
  # ThetaTrue: (p, q)
  p <- dim(ThetaEst)[1]
  q <- dim(ThetaEst)[2]
  L <- length(ThetaEst) / p / q
  ThetaEst <- array(ThetaEst, dim=c(p,q,L))
  rmse1 <- sapply(1:q, \(k) {
    diffTheta <- sweep(ThetaEst[, k, ], 1, phiTrueFunc$coefs[, k], "-")
    diffPhi <- B %*% diffTheta
    sqrt(colMeans(diffPhi^2))
  })
  rmse2 <- sapply(1:q, \(k) {
    diffTheta <- sweep(ThetaEst[, k, ], 1, phiTrueFunc$coefs[, k], "+")
    diffPhi <- B %*% diffTheta
    sqrt(colMeans(diffPhi^2))
  })
  pmin(rmse1, rmse2)
}
