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
  stepsize = 1e-3,
  maxIter = 1000,
  conv.check = FALSE,
  nIter.constStepSize = 0,
  nIter.adam = floor(maxIter),
  stepsize.decayrate = 0.5,
  stepsize.min = 0,
  nIter.slowerdecay = floor(0.5 * maxIter),
  stepsize.decayrate.slow = 0,
  beta1 = 0.9,
  beta2 = 0.99,
  adam.correction = TRUE,
  asgd.use = FALSE,
  asgd.start = 1,
  humblestart = 10,
  coord.scaling = FALSE,
  coord.scaling.start = 1,
  scales.cap.factor = 5,
  scales.cap.end = 200,
  nIter.1stTune = floor(0.1 * maxIter),
  nIter.lastTune = maxIter,
  period.tune = 100,
  vcrit.tune = c("ewmabv", "av"),
  ewmabv.beta = NULL,
  abv_aofv_w = 0.5,
  nIter.tauNoIncrease = floor(0.6 * maxIter),
  period.message = 100,
  verbose = TRUE,
  recordParams = TRUE,
  period.record = 20,
  period.time = period.record,
  test.mode = FALSE
) {
  ## Inputs:
  # TODO: obsGrid should be an attribute of data_generator
  # nIter.constStepSize: Number of steps with constant step size. Diminishing step sizes after that
  # nIter.adam: Number of Adam steps. Regular SGD after that

  ## Set params and initializations
  inits <- setParams.inits(inits, meanfun)
  optns.tau <- setParams.tau(tau, tau.control)
  tau <- sort(optns.tau$tau, decreasing = TRUE)
  taurange0 <- ifelse(length(tau) > 1, diff(range(log(tau))), 1)
  tau.control <- optns.tau$tau.control
  ntau <- tau.control$ntau
  adatau <- tau.control$adatau
  nselect <- tau.control$nselect
  nchild <- allocate_nchild(ntau, nselect)
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
  if (asgd.use) {
    vcrit.history$bv.avg <- array(dim = c(ntau, ntune))
    vcrit.history$av.avg <- array(dim = c(ntau, ntune))
    vcrit.history$ewmabv.avg <- array(dim = c(ntau, ntune))
  }
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
    theta_append = manifold.Stiefel.retractSinglePC(
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
    if (asgd.use) {
      params.history$Theta.avg <- params.history$Theta
      params.history$lambda.avg <- params.history$lambda
      params.history$sigma2.avg <- params.history$sigma2
    }
  } else {
    params.history <- NULL
  }
  iter.params <- 1

  # TODO: implement online mean estimation
  if (!is.null(theta_mu)) {
    stop("NOT IMPLEMENTED")
  }

  # stored stats for adaptive stochastic gradient
  if (all(c("grad_Theta", "grad_eta", "grad_zeta") %in% names(inits))) {
    inits$grad_Theta <- as.matrix(manifold.Stiefel.project(
      inits$grad_Theta,
      inits$Theta,
      G
    ))
    grad_all <- c(inits$grad_Theta, inits$grad_eta, inits$grad_zeta)
    m <- matrix(grad_all, nrow = p * q + q + 1, ncol = ntau)
    v <- matrix(grad_all^2, nrow = p * q + q + 1, ncol = ntau)
    # TODO: not finished. consider the case with initial grads
    # v2 <- c(
    #   colSums(inits$grad_Theta * G %*% inits$grad_Theta),
    #   inits$grad_eta^2, inits$grad_zeta^2
    # )
    # covg <- list(
    #   phi = as.matrix(t(inits$grad_Theta) %*% G %*% inits$grad_Theta),
    #   other = outer(c(inits$grad_eta, inits$grad_zeta), c(inits$grad_eta, inits$grad_zeta))
    # )
  } else {
    m <- matrix(0, nrow = p * q + q + 1, ncol = ntau)
    v <- matrix(0, nrow = p * q + q + 1, ncol = ntau)
    v2 <- matrix(0, nrow = q + q + 1, ncol = ntau)
    momt <- matrix(0, nrow = p * q + q + 1, ncol = ntau)
    covg <- list(
      phi = array(0, dim = c(q, q, ntau)),
      other = array(0, dim = c(q + 1, q + 1, ntau))
    )
    covg2 <- list(
      phi = array(0, dim = c(q, q, ntau)),
      other = array(0, dim = c(q + 1, q + 1, ntau))
    )
  }
  av <- numeric(ntau) # averaged one-function validation score
  ewmabv <- numeric(ntau) # averaged block validation score
  curr_bv <- numeric(ntau) # current block validation score
  if (asgd.use) {
    Theta.avg <- array(0, dim = c(p, q, ntau))
    lambda.avg <- array(0, dim = c(q, ntau))
    sigma2.avg <- array(0, dim = c(ntau))
    av.avg <- numeric(ntau)
    ewmabv.avg <- numeric(ntau)
    curr_bv.avg <- numeric(ntau)
  }

  total_count <- 0
  beta1_2i = beta1
  beta2_2i = beta2

  conv_flag <- rep(FALSE, ntune)
  for (i in 1:maxIter) {
    if (verbose && (i %% period.message == 0)) {
      message("Iter ", i)
    }
    if (i > nIter.constStepSize) {
      if (i <= nIter.slowerdecay) {
        if (stepsize.decayrate.slow == 0) {
          decay <- 1 / (1 + log(i - nIter.constStepSize))
        } else {
          decay <- 1 / (i^stepsize.decayrate.slow)
        }
      } else {
        if (i == nIter.slowerdecay + 1) {
          iter0 <- ceiling(1 / decay^(1 / stepsize.decayrate))
        }
        decay <- 1 / (iter0 + i - nIter.slowerdecay)^stepsize.decayrate
      }
    } else {
      decay <- 1
    }

    if (asgd.use && i == asgd.start) {
      asgd.count <- 1
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
      if ((!asgd.use) || i < asgd.start) {
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
      if ((!asgd.use) || i < asgd.start) {
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
        if ((!asgd.use) || i <= asgd.start) {
          vscore <- eval(parse(text = vcrit.tune))
        } else {
          vscore <- eval(parse(text = paste0(vcrit.tune, ".avg")))
        }
        lstar <- order(vscore)[1:nselect]
        parentId <- rep(lstar, nchild)
        tau.selectId <- unname(cbind(tau.selectId, parentId))
        m <- m[, parentId, drop = F]
        v <- v[, parentId, drop = F]
        v2 <- v2[, parentId, drop = F]
        momt <- momt[, parentId, drop = F]
        covg[['phi']] <- covg[['phi']][, , parentId, drop = F]
        covg[['other']] <- covg[['other']][, , parentId, drop = F]
        covg2[['phi']] <- covg2[['phi']][, , parentId, drop = F]
        covg2[['other']] <- covg2[['other']][, , parentId, drop = F]
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
        if (asgd.use && i >= asgd.start) {
          Theta.avg <- Theta.avg[,, parentId, drop = F]
          lambda.avg <- lambda.avg[, parentId, drop = F]
          sigma2.avg <- sigma2.avg[parentId]
          av.avg <- av.avg[parentId]
          ewmabv.avg <- ewmabv.avg[parentId]
        }
      }
    }

    if (asgd.use && i >= asgd.start) {
      asgd.count <- asgd.count + 1
    }

    for (l in 1:ntau) {
      time.start <- Sys.time()

      # if (i > nIter.lastTune && l > 1) break
      objective <- objfun(
        Ly,
        Ltid,
        asl(Theta, l),
        lambda[, l],
        sigma2[l],
        theta_mu,
        tau[l],
        stats = "grad"
      )
      grad_Theta <- objective$grad_Theta
      grad_Theta <- as.matrix(manifold.Stiefel.project(
        grad_Theta,
        asl(Theta, l),
        G
      ))
      grad_eta <- objective$grad_eta
      grad_zeta <- objective$grad_zeta
      grad_all <- c(grad_Theta, grad_eta, grad_zeta)
      if (verbose && l == 1 && (i %% period.message == 0)) {
        grad_norms <- sqrt(c(
          Theta = colSums(grad_Theta * G %*% grad_Theta),
          eta = sum(grad_eta^2), zeta = grad_zeta^2
        ))
        print(grad_norms)
      }

      m[, l] <- m[, l] * beta1 + grad_all * (1 - beta1)
      momt[1:(p*q), l] <- if (i == 1) {
        grad_Theta
      } else {
        beta1 * matrix(momt[1:(p*q), l], p, q) |>
          manifold.Stiefel.project(asl(Theta, l), G) +
          (1 - beta1) * grad_Theta
      }
      momt[(p*q+1):(p*q+q+1), l] <- momt[(p*q+1):(p*q+q+1), l] * beta1 +
        (1 - beta1) * c(grad_eta, grad_zeta)
      v[, l] <- v[, l] * beta2 + grad_all^2 * (1 - beta2)
      v2[, l] <- v2[, l] * (i - 1) / i + 1 / i * c(
        colSums(grad_Theta * G %*% grad_Theta),
        grad_eta^2, grad_zeta^2
      )
      covg[['phi']][,,l] <- covg[['phi']][,,l] * (i - 1) / i +
        1 / i * as.matrix(t(grad_Theta) %*% G %*% grad_Theta)
      covg[['other']][,,l] <- covg[['other']][,,l] * (i - 1) / i +
        1 / i * outer(c(grad_eta, grad_zeta), c(grad_eta, grad_zeta))
      covg2[['phi']][,,l] <- covg2[['phi']][,,l] * (i - 1) / i +
        1 / i * as.matrix(t(grad_Theta) %*% G %*% grad_Theta)
      covg2[['other']][,,l] <- covg2[['other']][,,l] * beta2 +
        (1 - beta2) * outer(c(grad_eta, grad_zeta), c(grad_eta, grad_zeta))
      if (i <= nIter.adam) {
        # adaptive gradients
        direc <- m[, l] / (sqrt(v[, l]) + 1e-8)
      } else {
        # switch back to regular sgd
        # TODO: modified the weighting rule
        # (1) usual normalized SGD
        # direc <- grad_all / (sqrt(mean(grad_all^2)) + 1e-8)
        # (2) an intuitive AdaGrad
        # direc <- c(
        #   sweep(grad_Theta, 2, sqrt(v2[1:q,l]), "/"),
        #   c(grad_eta, grad_zeta) / sqrt(tail(v2[,l], q+1))
        # )
        # (3) another intuitive AdaGrad
        nugget <- ifelse(i < q+1, 1, 0)
        covg_sqinv_phi <- pracma::sqrtm(covg2[['phi']][,,l] + diag(1e-6+nugget,q,q))$Binv
        covg_sqinv_other <- pracma::sqrtm(covg2[['other']][,,l] + diag(1e-6+nugget,q+1,q+1))$Binv
        direc <- c(
          grad_Theta %*% covg_sqinv_phi,
          covg_sqinv_other %*% c(grad_eta, grad_zeta)
        )
        # (4) an intuitive Adam
        # nugget <- ifelse(i < q+1, 1, 0)
        # covg_sqinv_phi <- pracma::sqrtm(covg2[['phi']][,,l] + diag(1e-6+nugget,q,q))$Binv
        # covg_sqinv_other <- pracma::sqrtm(covg2[['other']][,,l] + diag(1e-6+nugget,q+1,q+1))$Binv
        # direc <- c(
        #   grad_Theta %*% covg_sqinv_phi,
        #   covg_sqinv_other %*% c(grad_eta, grad_zeta)
        # )
        if (adam.correction) {
          direc <- direc * (sqrt(1 - beta2_2i) / (1 - beta1_2i))
        }
      }

      if (verbose && l == 1 && (i %% period.message == 0)) {
        direc_norms <- sqrt(c(
          Theta = colSums(matrix(direc[1:(p*q)],p,q) * G %*% matrix(direc[1:(p*q)],p,q)),
          eta = sum(direc[(p*q+1):(p*q+q)]^2),
          zeta = tail(direc,1)^2
        ))
        print(direc_norms)
      }

      if (coord.scaling && i > coord.scaling.start) {
        scales.cap <- (1 / scales.cap.factor^(0:(q - 1)))
        if (i <= scales.cap.end) {
          if (asgd.use && i >= asgd.start) {
            lambda.scales <- pmax(lambda.avg, lambda.avg[1] * scales.cap)
          } else {
            lambda.scales <- pmax(lambda, lambda[1] * scales.cap)
          }
        } else {
          lambda.scales <- lambda
        }
        if (asgd.use && i >= asgd.start) {
          coord.scales <- 1 /
            c(rep(2 * lambda.scales, each = p), lambda.scales, sigma2.avg)
        } else {
          coord.scales <- 1 /
            c(rep(2 * lambda.scales, each = p), lambda.scales, sigma2)
        }
        coord.scales <- coord.scales / coord.scales[1]
        direc <- direc * coord.scales
      }

      # update
      humble_factor <- ifelse(i <= humblestart, 0.1, 1)
      curr_stepsize <- humble_factor * stepsize * decay
      if (curr_stepsize < stepsize.min) {
        curr_stepsize <- stepsize.min
      }
      Theta[,, l] <- Theta[,, l] - curr_stepsize * direc[1:(p * q)]
      lambda[, l] <- lambda[, l] *
        exp(-curr_stepsize * direc[(p * q + 1):(p * q + q)])
      sigma2[l] <- sigma2[l] * exp(-curr_stepsize * direc[p * q + q + 1])

      if (
        all(
          curr_stepsize * abs(direc) <
            1e-8 * abs(c(Theta[,, l], log(lambda[, l]), log(sigma2[l])))
        )
      ) {
        conv_flag[l] <- TRUE
      }

      if (asgd.use && i >= asgd.start) {
        Theta.avg[,, l] <- Theta.avg[,, l] *
          ((asgd.count - 1) / asgd.count) +
          Theta[,, l] / asgd.count
        lambda.avg[, l] <- lambda.avg[, l]^((asgd.count - 1) / asgd.count) *
          lambda[, l]^(1 / asgd.count)
        sigma2.avg[l] <- sigma2.avg[l]^((asgd.count - 1) / asgd.count) *
          sigma2[l]^(1 / asgd.count)
      }

      # retraction
      Theta[,, l] <- manifold.Stiefel.retract(asl(Theta, l), G)
      if (asgd.use) {
        if (i >= asgd.start) {
          Theta.avg[,, l] <- manifold.Stiefel.retract(asl(Theta.avg, l), G)
        } else {
          Theta.avg[,, l] <- Theta[,, l]
          lambda.avg[, l] <- lambda[, l]
          sigma2.avg[l] <- sigma2[l]
        }
      }

      # parallel transport of the adam m
      m[1:(p * q), l] <- as.matrix(
        manifold.Stiefel.project(
          matrix(m[1:(p * q), l], p, q),
          asl(Theta, l),
          G
        )
      )

      # record time
      time.end <- Sys.time()
      time.history[l, i] <- time.history[l, i] +
        difftime(time.end, time.start, units = 'secs')
    }

    # record parameters
    if (recordParams && (i %% period.record == 0)) {
      iter.params <- c(iter.params, i)
      idx <- i %/% period.record + 1
      params.history$Theta[,,, idx] <- Theta
      params.history$lambda[,, idx] <- lambda
      params.history$sigma2[, idx] <- sigma2
      if (asgd.use) {
        params.history$Theta.avg[,,, idx] <- Theta.avg
        params.history$lambda.avg[,, idx] <- lambda.avg
        params.history$sigma2.avg[, idx] <- sigma2.avg
      }
    }

    beta1_2i = beta1_2i * beta1
    beta2_2i = beta2_2i * beta2

    if (conv.check && all(conv_flag)) break
  }

  # Final selection
  if (!asgd.use) {
    minId <- which.min(ewmabv)
    tau.min <- tau[minId]
  } else {
    minId <- which.min(ewmabv.avg)
    tau.min <- tau[minId]
  }

  out <- list(
    Theta = Theta,
    lambda = lambda,
    sigma2 = sigma2,
    theta_mu = theta_mu,
    tau = tau,
    tau.min = tau.min,
    av = av,
    ewmabv = ewmabv,
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
    asgd.use = asgd.use
  )

  if (asgd.use) {
    out <- append(
      out,
      list(
        Theta.avg = Theta.avg,
        lambda.avg = lambda.avg,
        sigma2.avg = sigma2.avg,
        av.avg = av.avg,
        ewmabv.avg = ewmabv.avg
      )
    )
  }

  return(out)
}

# fpca.rasa <- function(
#   data_generator,
#   obsGrid,
#   inits,
#   meanfun = FALSE,
#   tau = NULL,
#   tau.control = list(),
#   nbatch = 1,
#   stepsize = 1e-3,
#   maxIter = 1000,
#   conv.check = FALSE,
#   stepsize.decayrate = 0.5,
#   stepsize.min = 0,
#   nIter.slowerdecay = floor(0.5 * maxIter),
#   stepsize.decayrate.slow = 0,
#   asgd.start = 1,
#   beta0 = 0.9,
#   beta1 = 0.9,
#   beta2 = 0.99,
#   nIter.1stTune = floor(0.1 * maxIter),
#   nIter.lastTune = maxIter,
#   period.tune = 100,
#   vcrit.tune = c("ewmabv", "av"),
#   ewmabv.beta = NULL,
#   abv_aofv_w = 0.5,
#   nIter.tauNoIncrease = floor(0.6 * maxIter),
#   period.message = 100,
#   verbose = TRUE,
#   recordParams = TRUE,
#   period.record = 20,
#   period.time = period.record,
#   test.mode = FALSE
# ) {
#   # FIXME: FUNCTION OUTDATED
#   ## Inputs:
#   # nIter.constStepSize: Number of steps with constant step size. Diminishing step sizes after that
#   # nIter.adam: Number of Adam steps. Regular SGD after that

#   ## Set params and initializations
#   inits <- setParams.inits(inits, meanfun)
#   optns.tau <- setParams.tau(tau, tau.control)
#   tau <- sort(optns.tau$tau, decreasing = TRUE)
#   taurange0 <- ifelse(length(tau) > 1, diff(range(log(tau))), 1)
#   tau.control <- optns.tau$tau.control
#   ntau <- tau.control$ntau
#   adatau <- tau.control$adatau
#   nselect <- tau.control$nselect
#   nchild <- allocate_nchild(ntau, nselect)
#   if (adatau) {
#     tau.history <- tau
#     tau.selectId <- c()
#   } else {
#     tau.history <- NULL
#     tau.selectId <- NULL
#   }
#   ntune <- floor((maxIter - nIter.1stTune) / period.tune) + 1
#   if (ntune < 5) {
#     ewmabv.beta <- 0.5
#   } else {
#     ewmabv.beta.candidates <- seq(0.1, 0.95, 0.05)
#     valid.ewmabv.beta.id <- which(
#       ewmabv.beta.candidates^(ntune - 1) < 0.5 / ntune^1.2
#     )
#     if (length(valid.ewmabv.beta.id) == 0) {
#       idx <- 1
#       warning(
#         "Improper choice of ewmabv.beta. The tuning period is possibly too large."
#       )
#     }
#     idx <- max(valid.ewmabv.beta.id)
#     if (is.null(ewmabv.beta) || ewmabv.beta > ewmabv.beta.candidates[idx]) {
#       ewmabv.beta <- ewmabv.beta.candidates[idx]
#     }
#   }
#   vcrit.tune <- match.arg(vcrit.tune)
#   vcrit.history <- list(
#     av = array(dim = c(ntau, ntune)),
#     ewmabv = array(dim = c(ntau, ntune))
#   )
#   vcrit.history$av.avg <- array(dim = c(ntau, ntune))
#   vcrit.history$ewmabv.avg <- array(dim = c(ntau, ntune))
#   iter.select <- 1
#   stopifnot(stepsize > 0)
#   stopifnot(stepsize.decayrate > 0 && stepsize.decayrate <= 1)
#   stopifnot(stepsize.decayrate.slow >= 0 && stepsize.decayrate.slow <= 1)

#   p <- nrow(inits$Theta)
#   q <- ncol(inits$Theta)
#   Theta <- array(inits$Theta, dim = c(p, q, ntau))
#   lambda <- array(inits$lambda, dim = c(q, ntau))
#   sigma2 <- array(inits$sigma2, dim = c(ntau))
#   theta_mu <- inits$theta_mu

#   if (recordParams) {
#     nrecord <- floor(maxIter / period.record) + 1
#     params.history <- list(
#       Theta = array(dim = c(p, q, ntau, nrecord)),
#       lambda = array(dim = c(q, ntau, nrecord)),
#       sigma2 = array(dim = c(ntau, nrecord))
#     )
#     params.history$Theta[,,, 1] <- Theta
#     params.history$lambda[,, 1] <- lambda
#     params.history$sigma2[, 1] <- sigma2
#     params.history$Theta.avg <- params.history$Theta
#     params.history$lambda.avg <- params.history$lambda
#     params.history$sigma2.avg <- params.history$sigma2
#   } else {
#     params.history <- NULL
#   }
#   iter.params <- 1

#   # TODO: implement online mean estimation
#   if (!is.null(theta_mu)) {
#     stop("NOT IMPLEMENTED")
#   }

#   # matrix part: RASA algorithm
#   m.Theta <- matrix(0, nrow = p * q, ncol = ntau)
#   rasa.l <- matrix(0, nrow = p, ncol = ntau)
#   rasa.r <- matrix(0, nrow = q, ncol = ntau)
#   # vector part: regular Adam
#   m.lamsig <- matrix(0, nrow = q + 1, ncol = ntau)
#   v.lamsig <- matrix(0, nrow = q + 1, ncol = ntau)

#   av <- numeric(ntau) # averaged one-function validation score
#   ewmabv <- numeric(ntau) # averaged block validation score
#   curr_bv <- numeric(ntau) # current block validation score
#   Theta.avg <- array(0, dim = c(p, q, ntau))
#   lambda.avg <- array(0, dim = c(q, ntau))
#   sigma2.avg <- array(0, dim = c(ntau))
#   av.avg <- numeric(ntau)
#   ewmabv.avg <- numeric(ntau)
#   curr_bv.avg <- numeric(ntau)

#   total_count <- 0
#   time <- c()
#   time.start <- Sys.time()

#   conv_flag <- rep(FALSE, ntune)
#   for (i in 1:maxIter) {
#     if (verbose && (i %% period.message == 0)) {
#       message("Iter ", i)
#     }
#     if (i <= nIter.slowerdecay) {
#       if (stepsize.decayrate.slow == 0) {
#         decay <- 1 / (1 + log(i))
#       } else {
#         decay <- 1 / (i^stepsize.decayrate.slow)
#       }
#     } else {
#       if (i == nIter.slowerdecay + 1) {
#         iter0 <- ceiling(1 / decay^(1 / stepsize.decayrate))
#       }
#       decay <- 1 / (iter0 + i - nIter.slowerdecay)^stepsize.decayrate
#     }

#     if (i == asgd.start) {
#       asgd.count <- 1
#       Theta.avg <- Theta
#       lambda.avg <- lambda
#       sigma2.avg <- sigma2
#     }

#     # receive new data
#     obs <- data_generator(nbatch, total_count)
#     Ltid <- obs$Ltid
#     Ly <- obs$Ly
#     total_count <- total_count + nbatch

#     # calculate validation score
#     for (l in 1:ntau) {
#       # if (i > nIter.lastTune && l > 1) break
#       if (i < asgd.start) {
#         lik <- objfun(
#           Ly,
#           Ltid,
#           asl(Theta, l),
#           lambda[, l],
#           sigma2[l],
#           theta_mu,
#           tau[l],
#           stats = "loss"
#         )$lik
#         av[l] <- av[l] * (i - 1) / i + lik * 1 / i
#         curr_bv[l] <- curr_bv[l] + lik
#       } else {
#         lik.avg <- objfun(
#           Ly,
#           Ltid,
#           asl(Theta.avg, l),
#           lambda.avg[, l],
#           sigma2.avg[l],
#           theta_mu,
#           tau[l],
#           stats = "loss"
#         )$lik
#         av.avg[l] <- av.avg[l] * (i - 1) / i + lik.avg * 1 / i
#         curr_bv.avg[l] <- curr_bv.avg[l] + lik.avg
#       }
#     }

#     if (
#       i > 1 && i >= nIter.1stTune && (i - nIter.1stTune) %% period.tune == 0
#     ) {
#       itune <- floor((i - nIter.1stTune) / period.tune) + 1
#       if (itune == 1) {
#         curr_bv <- curr_bv / (i * nbatch)
#       } else {
#         curr_bv <- curr_bv / (period.tune * nbatch)
#       }
#       if (i < asgd.start) {
#         ewmabv <- ewmabv * ewmabv.beta + curr_bv * (1 - ewmabv.beta)
#         curr_bv <- numeric(ntau)
#         vcrit.history$av[, itune] <- av
#         vcrit.history$ewmabv[, itune] <- ewmabv
#       } else {
#         curr_bv.avg <- curr_bv.avg / (period.tune * nbatch)
#         ewmabv.avg <- ewmabv.avg * ewmabv.beta + curr_bv.avg * (1 - ewmabv.beta)
#         curr_bv.avg <- numeric(ntau)
#         vcrit.history$av.avg[, itune] <- av.avg
#         vcrit.history$ewmabv.avg[, itune] <- ewmabv.avg
#       }

#       # online selection of tau
#       if (adatau && i <= nIter.lastTune) {
#         iter.select <- c(iter.select, i)
#         if (i <= asgd.start) {
#           vscore <- eval(parse(text = vcrit.tune))
#         } else {
#           vscore <- eval(parse(text = paste0(vcrit.tune, ".avg")))
#         }
#         lstar <- order(vscore)[1:nselect]
#         parentId <- rep(lstar, nchild)
#         tau.selectId <- unname(cbind(tau.selectId, parentId))
#         rasa.l <- rasa.l[, parentId, drop = F]
#         rasa.r <- rasa.r[, parentId, drop = F]
#         m <- m[, parentId, drop = F]
#         v <- v[, parentId, drop = F]
#         av <- av[parentId]
#         ewmabv <- ewmabv[parentId]
#         Theta <- Theta[,, parentId, drop = F]
#         lambda <- lambda[, parentId, drop = F]
#         sigma2 <- sigma2[parentId]
#         tau.control <- setParams.tau(
#           tau,
#           tau.control,
#           delta.min = taurange0 / max(1, ntau - 2) / (itune^0.7)
#         )$tau.control
#         tau <- explore_newtau(
#           tau[lstar],
#           nchild,
#           tau.control$delta,
#           i > nIter.tauNoIncrease
#         )
#         tau.history <- unname(cbind(tau.history, tau))
#         if (i >= asgd.start) {
#           Theta.avg <- Theta.avg[,, parentId, drop = F]
#           lambda.avg <- lambda.avg[, parentId, drop = F]
#           sigma2.avg <- sigma2.avg[parentId]
#           av.avg <- av.avg[parentId]
#           ewmabv.avg <- ewmabv.avg[parentId]
#         }
#       }
#     }

#     if (i >= asgd.start) {
#       asgd.count <- asgd.count + 1
#     }

#     for (l in 1:ntau) {
#       # if (i > nIter.lastTune && l > 1) break
#       objective <- objfun(
#         Ly,
#         Ltid,
#         asl(Theta, l),
#         lambda[, l],
#         sigma2[l],
#         theta_mu,
#         tau[l],
#         stats = "grad"
#       )
#       grad_Theta <- objective$grad_Theta
#       grad_Theta <- as.matrix(manifold.Stiefel.project(
#         grad_Theta,
#         asl(Theta, l),
#         G
#       ))
#       grad_eta <- objective$grad_eta
#       grad_zeta <- objective$grad_zeta
#       grad_all <- c(grad_Theta, grad_eta, grad_zeta)

#       m.Theta[, l] <- m.Theta[, l] * beta1 + grad_Theta * (1 - beta1)
#       rasa.l[, l] <- rasa.l[, l] *
#         beta2 +
#         c(diag(grad_Theta %*% t(grad_Theta))) * (1 - beta2) / q
#       rasa.r[, l] <- rasa.r[, l] *
#         beta2 +
#         c(diag(t(grad_Theta) %*% grad_Theta)) * (1 - beta2) / p
#       m.lamsig[, l] <- m.lamsig[, l] *
#         beta1 +
#         c(grad_eta, grad_zeta) * (1 - beta1)
#       v.lamsig[, l] <- v.lamsig[, l] *
#         beta2 +
#         c(grad_eta, grad_zeta)^2 * (1 - beta2)
#       direc.Theta <- outer(
#         (rasa.l[, l] + 1e-8)^(-1 / 4),
#         (rasa.r[, l] + 1e-8)^(-1 / 4)
#       ) *
#         matrix(m.Theta[, l], p, q)
#       direc.Theta <- as.matrix(manifold.Stiefel.project(
#         direc.Theta,
#         asl(Theta, l),
#         G
#       ))
#       direc.lamsig <- m.lamsig[, l] / (sqrt(v.lamsig[, l]) + 1e-8)
#       direc <- c(direc.Theta, direc.lamsig)

#       # update
#       curr_stepsize <- stepsize * decay
#       if (curr_stepsize < stepsize.min) {
#         curr_stepsize <- stepsize.min
#       }
#       Theta[,, l] <- Theta[,, l] - curr_stepsize * direc[1:(p * q)]
#       lambda[, l] <- lambda[, l] *
#         exp(-curr_stepsize * direc[(p * q + 1):(p * q + q)])
#       sigma2[l] <- sigma2[l] * exp(-curr_stepsize * direc[p * q + q + 1])

#       if (
#         all(
#           curr_stepsize * abs(direc) <
#             1e-8 * abs(c(Theta[,, l], log(lambda[, l]), log(sigma2[l])))
#         )
#       ) {
#         conv_flag[l] <- TRUE
#       }

#       if (i >= asgd.start) {
#         Theta.avg[,, l] <- Theta.avg[,, l] *
#           ((asgd.count - 1) / asgd.count) +
#           Theta[,, l] / asgd.count
#         lambda.avg[, l] <- lambda.avg[, l]^((asgd.count - 1) / asgd.count) *
#           lambda[, l]^(1 / asgd.count)
#         sigma2.avg[l] <- sigma2.avg[l]^((asgd.count - 1) / asgd.count) *
#           sigma2[l]^(1 / asgd.count)
#       }

#       # retraction
#       Theta[,, l] <- manifold.Stiefel.retract(asl(Theta, l), G)
#       if (i >= asgd.start) {
#         Theta.avg[,, l] <- manifold.Stiefel.retract(asl(Theta.avg, l), G)
#       } else {
#         Theta.avg[,, l] <- Theta[,, l]
#         lambda.avg[, l] <- lambda[, l]
#         sigma2.avg[l] <- sigma2[l]
#       }
#     }

#     if (i %% period.time == 0) {
#       time.end <- Sys.time()
#       time <- c(time, difftime(time.end, time.start, units = 'secs'))
#       time.start <- Sys.time()
#     }

#     # record parameters
#     if (recordParams && (i %% period.record == 0)) {
#       iter.params <- c(iter.params, i)
#       idx <- i %/% period.record + 1
#       params.history$Theta[,,, idx] <- Theta
#       params.history$lambda[,, idx] <- lambda
#       params.history$sigma2[, idx] <- sigma2
#       params.history$Theta.avg[,,, idx] <- Theta.avg
#       params.history$lambda.avg[,, idx] <- lambda.avg
#       params.history$sigma2.avg[, idx] <- sigma2.avg
#     }

#     if (conv.check && all(conv_flag)) break
#   }

#   # Final selection
#   minId <- which.min(ewmabv.avg)
#   tau.min <- tau[minId]

#   out <- list(
#     Theta = Theta,
#     lambda = lambda,
#     sigma2 = sigma2,
#     theta_mu = theta_mu,
#     tau = tau,
#     tau.min = tau.min,
#     av = av,
#     ewmabv = ewmabv,
#     tau.select = list(
#       tau.history = tau.history,
#       tau.selectId = tau.selectId,
#       iter.select = iter.select
#     ),
#     params.history = list(params = params.history, iter.params = iter.params),
#     vcrit.history = vcrit.history,
#     time = time,
#     stepsize = stepsize,
#     stepsize.final = stepsize * decay,
#     ewmabv.beta = ewmabv.beta,
#     Theta.avg = Theta.avg,
#     lambda.avg = lambda.avg,
#     sigma2.avg = sigma2.avg,
#     av.avg = av.avg,
#     ewmabv.avg = ewmabv.avg
#   )

#   return(out)
# }

objfun <- function(
  Ly,
  Ltid,
  Theta,
  lambda,
  sigma2,
  theta_mu = NULL,
  tau = 1e-3,
  stats = c("all", "loss", "grad")
) {
  # global constants: B, Omega
  stats <- match.arg(stats)
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
        1 / n * (
          1 / sigma2 *
            (sum(Yi^2) - sum(PhiT_Yi * invQ %*% PhiT_Yi)) +
            (mi - q) * log(sigma2) +
            sum(log(lambda)) +
            logdetQ
        )
    }
    if (stats == "all" || stats == "grad") {
      Phit_Phi <- crossprod(Phi)
      invQ_PhiT_Yi <- invQ %*% PhiT_Yi
      grad_Theta <- grad_Theta +
        2 / n *
          as.matrix(
            crossprod(Bi, Phi) %*% (
              1 / sigma2 * invQ_PhiT_Yi %*% t(invQ_PhiT_Yi)
                + invQ
            )
              - 1 / sigma2 * crossprod(Bi, Yi) %*% t(invQ_PhiT_Yi)
          )

      grad_eta <- grad_eta +
        1 / n * (
          -invQ_PhiT_Yi^2 / lambda + rep(1, q)
            - sigma2 * diag(invQ) / lambda
        )

      grad_zeta <- grad_zeta +
        1 / n * (
          -sum(Yi^2) / sigma2
            + sum(invQ_PhiT_Yi * PhiT_Yi) / sigma2
            + sum(invQ_PhiT_Yi^2 / lambda)
            + mi - q
            + sigma2 * sum(diag(invQ) / lambda)
        )
    }
  }
  out <- list()
  if (stats == "all" || stats == "loss") {
    out[["fval"]] <- fval + sum((Theta * Omega %*% Theta) %*% tau)
    out[["lik"]] <- fval # unpenalized likelihood
  }
  if (stats == "all" || stats == "grad") {
    # TODO: calculation changed. Check performance
    out[["grad_Theta"]] <- solve(
      G, grad_Theta +
        2 * matrix(rep(tau, each = p), nrow = p) * as.matrix(Omega %*% Theta)
    ) |> as.matrix()
    # out[["grad_Theta"]] <- grad_Theta +
    #     2 * matrix(rep(tau, each = p), nrow = p) * as.matrix(Omega %*% Theta)
    out[["grad_eta"]] <- as.numeric(grad_eta)
    out[["grad_zeta"]] <- as.numeric(grad_zeta)
  }
  return(out)
}



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

setParams.inits <- function(inits, meanfun = FALSE) {
  if (!("Theta" %in% names(inits))) {
    stop("Lack an initialization for Theta")
  }
  inits$Theta <- manifold.Stiefel.retract(inits$Theta, G)
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

setParams.tau <- function(tau = NULL, tau.control = list(), delta.min = NULL) {
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
    if (is.null(delta.min)) {
      delta.min <- 0.3
    } else {
      delta.min <- max(delta.min, 0.3)
    }
    if (length(tau) > 1) {
      tau.control$delta <- max(delta.min, mean(diff(sort(log(tau)))) / 3)
    } else {
      tau.control$delta <- 0
    }
  }
  return(list(tau = tau, tau.control = tau.control))
}


# setParams.sgd <- function(sgd.control) {
#   # stepsize = 1e-3,
#   # maxIter = 1000,
#   # nIter.constStepSize = 0, nIter.adam = floor(maxIter),
#   # stepsize.decayrate = 0.5, stepsize.min = 0,
#   # nIter.slowerdecay = floor(0.5 * maxIter),
#   # stepsize.decayrate.slow = 0,
#   # beta1 = 0.9, beta2 = 0.99,
#   # asgd.use = FALSE, asgd.start = 100
#   if (! "stepsize" %in% sgd.control) stepsize <- 1e-3
#   if (! "maxIter" %in% sgd.control) maxIter <- 1000
#   if (! "nIter.constStepSize" %in% sgd.control) nIter.constStepSize <- 0
# }

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
