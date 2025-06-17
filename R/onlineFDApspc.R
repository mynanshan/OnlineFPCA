fpca.pspc.online <- function(
    data_generator, accumulate_data = FALSE,
    basis = NULL, init = NULL,
    nbasis = 10, npc = NULL,
    tau_mu = 10^(-4), tau = 10^(0),
    dlogtau = 1, nexplore = 3, dlogtau_decay = 0.2,
    n = 20, step_size = 1e-3, momentum = 0.9,
    numRound = 10, numIterPerRound = 40, maxRoundLinesearch = NULL,
    maxIterLinesearch = 20, minIterLinesearch = 0,
    step_size_enlarge_factor = 0.1,
    armijo.beta = 0.5, armijo.delta = 0.4,
    rasa.beta = 0.99,
    record.estimate = FALSE,
    record.runtime = FALSE,
    stream_seeds = NULL,
    verbose = TRUE) {
  
  # linesearch_rasa: rasa after linesearch
  # vanilla: vanilla RSGD
  # linesearch: RSGD with ad-hoc line search
  # rasa: Riemannian Adam by (Kasai, 2019)
  
  if (is.null(basis)) {
    stop("NOT IMPLEMENTED")
    init <- NULL
  }
  stopifnot(class(basis) %in% c("basisfd", "TensorBasis"))
  if (is(basis, "basisfd")) {
    p <- basis$nbasis
  } else if (is(basis, "TensorBasis")) {
    p <- attr(basis, "nbasis")
  }
  
  if (is.null(init)) stop("NOT IMPLEMENTED")
  if (!is.null(init)) {
    stopifnot(all(c("Theta", "theta_mu", "sigma2", "lambda") %in% names(init)))
    Theta <- init[['Theta']]
    theta_mu <- init[['theta_mu']]
    sigma2 <- init[['sigma2']]
    lambda <- init[['lambda']]
    mean_Bt_B <- init[['mean_Bt_B']]
    mean_Bt_Y <- init[['mean_Bt_Y']]
    N <- init[['N']]
  }
  
  
  q <- npc
  if (length(lambda) < q) {
    lambda <- c(lambda, rep(min(lambda) * 0.1, q - length(lambda)))
  } else lambda <- lambda[1:q]
  
  G <- get_basis_inprod_matrix(basis)
  Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
  
  if (numIterPerRound < maxIterLinesearch) {
    warning("Reset numIterPerRound to maxIterLinesearch")
    numIterPerRound <- maxIterLinesearch
  }
  
  if (record.estimate) {
    est_trace <- list(
      Theta = c(Theta), theta_mu = theta_mu,
      lambda = lambda, sigma2 = sigma2, avg_mse = c()
    )
  }
  
  if (nexplore %% 2 == 0) nexplore <- nexplore + 1
  tau_parent <- tau
  ntau_parent <- length(tau)
  ndexpl <- seq((nexplore - 1) / 2, -(nexplore - 1) / 2, -1)
  ntau <- ntau_parent * nexplore
  # replicates
  step_size <- rep(step_size, ntau_parent)
  step_size_enlarge_factor <- rep(step_size_enlarge_factor, ntau_parent)
  Theta_parent <- array(Theta, dim = c(p, q, ntau_parent))
  lambda_parent <- array(lambda, dim = c(q, ntau_parent))
  sigma2_parent <- rep(sigma2, ntau_parent)
  # to store the candidates
  Theta <- array(0, dim = c(p, q, ntau))
  lambda <- array(0, dim = c(q, ntau))
  sigma2 <- numeric(ntau)
  
  # maintained statistics
  if (is.null(init)) {
    mean_Bt_B <- matrix(0, p, p)
    mean_Bt_Y <- matrix(0, p, 1)
    N <- 0
    m_mean <- 0
    m2_mean <- 0
  }
  MSE <- numeric(ntau_parent)
  rasa.l <- rasa.lhat <- matrix(0, nrow = p, ncol = ntau_parent)
  rasa.r <- rasa.rhat <- matrix(0, nrow = q, ncol = ntau_parent)
  
  # traces of smoothness parameters
  tau_path <- list(list(tau = tau_parent))
  stepsize_path <- c(step_size)
  
  if (record.runtime) runtime <- c()
  if (accumulate_data) {
    cumLt <- list()
    cumLy <- list()
    cumLmi <- c()
  }
  
  # SGD updates ----------------
  niter <- 1
  for (ir in 1:numRound) {
    
    if (verbose) message("> Round ", ir)
    
    # determine numbers of iterations for Stage 1,2
    if (!is.null(maxRoundLinesearch) && ir > maxRoundLinesearch) {
      numIterStage1 <- 0
    } else {
      numIterStage1 <- max(round(maxIterLinesearch / ir), minIterLinesearch)
    }
    numIterStage2 <- numIterPerRound - numIterStage1
    
    ## initialize new candidates
    tau <- c(sapply(tau_parent,
                    \(x) x * exp(dlogtau / (ir^dlogtau_decay) * ndexpl)))
    for (ii in 1:ntau_parent) {
      tmpidx <- ((ii-1)*nexplore+1):(ii*nexplore)
      Theta[,,tmpidx] <- asl(Theta_parent,ii)
      sigma2[tmpidx] <- sigma2_parent[ii]
      lambda[,tmpidx] <- lambda_parent[,ii]
    }
    step_size <- rep(step_size, each = nexplore)
    step_size_enlarge_factor <- rep(step_size_enlarge_factor, each = nexplore)
    MSE <- rep(MSE, each = nexplore)
    rasa.l <- do.call(
      cbind, lapply(1:ntau_parent, \(ii) replicate(nexplore, rasa.l[,ii])))
    rasa.r <- do.call(
      cbind, lapply(1:ntau_parent, \(ii) replicate(nexplore, rasa.r[,ii])))
    rasa.lhat <- do.call(
      cbind, lapply(1:ntau_parent, \(ii) replicate(nexplore, rasa.lhat[,ii])))
    rasa.rhat <- do.call(
      cbind, lapply(1:ntau_parent, \(ii) replicate(nexplore, rasa.rhat[,ii])))
    
    ## main iterations: stochastic optimization
    
    for (it in 1:numIterPerRound) {
      
      # message(">> Iter ", it)
      
      if (!is.null(stream_seeds))
        set.seed(stream_seeds[niter %% length(stream_seeds)])
      obs <- data_generator(n)
      Lt <- obs$Lt; Ly <- obs$Ly; Lmi <- obs$Lmi
      m <- sum(Lmi)
      m2 <- sum(Lmi * (Lmi - 1))
      consts <- list(p = p, q = q, n = n, m2 = m2, Lmi = Lmi)
      
      time.start <- Sys.time()
      
      ## update mean ---------- 
      # TODO: add automatic tuning for mean function
      for (i in 1:n) {
        Yi <- Ly[[i]]; Ti <- Lt[[i]]; mi <- Lmi[i]
        Bi <- eval_basis(Ti, basis)
        mean_Bt_B <- mean_Bt_B * N / (N + n) +
          Matrix::crossprod(Bi) / (N + n)
        mean_Bt_Y <- mean_Bt_Y * N / (N + n) +
          Matrix::crossprod(Bi, Yi) / (N + n)
      }
      theta_mu <- as.matrix(Matrix::solve(mean_Bt_B + tau_mu * Omega, mean_Bt_Y))
      
      # gradients and MSE
      if (!accumulate_data) {
        optObj <- pspc_stats(
          Ly, Lt, basis, Theta, theta_mu, sigma2, lambda,
          stats = c("grad", "mse"))
      } else {
        cumLy <- c(cumLy, Ly)
        cumLt <- c(cumLt, Lt)
        optObj <- pspc_stats(
          cumLy, cumLt, basis, Theta, theta_mu, sigma2, lambda,
          stats = c("grad", "mse"))
      }
      
      grad_Theta <- sapply(1:ntau, \(itau) {
        asl(optObj$grad_Theta,itau) +
          2 * tau[itau] * as.matrix(Omega %*% asl(Theta,itau))
      }, simplify = "array")
      grad_Theta <- sapply(1:ntau, \(itau) {
        as.matrix(manifold.Stiefel.project(asl(grad_Theta,itau), asl(Theta,itau), G))
      }, simplify = "array")
      grad_zeta <- optObj$grad_zeta
      grad_eta <- optObj$grad_eta
      mse_curr <- optObj$mse
      ## update parameters
      step_decay <- 1 / sqrt(niter)
      # print(step_decay)
      if (it == 1) armijo.alpha <- step_size
      if (it <= numIterStage1) {
        # Stage 1: Adapt step sizes with line-search ----------------
        for (itau in 1:ntau) {
          # message(">>> Line Search ", itau)
          sg_obj <- sg.linesearch(
            asl(Theta,itau), theta_mu, sigma2[itau], lambda[,itau], tau[itau],
            asl(grad_Theta,itau), grad_zeta[itau], grad_eta[,itau],
            step_size_enlarge_factor[itau], step_decay,
            armijo.alpha[itau], armijo.beta, armijo.delta,
            Ly, Lt, basis, consts = c(consts, list(
              G = G, Omega = Omega, mse_curr = mse_curr[itau])
          )) 
          Theta[,,itau] <- sg_obj$Theta 
          sigma2[itau] <- sg_obj$sigma2
          lambda[,itau] <- sg_obj$lambda
          armijo.alpha[itau] <- sg_obj$armijo.alpha
          step_size_enlarge_factor[itau] <- sg_obj$step_size_enlarge_factor
          
          # maintain the stats for RASA in the meantime
          rasa.l[,itau] <- rasa.beta * rasa.l[,itau] +
            (1 - rasa.beta) * rowSums(asl(grad_Theta,itau)^2) / q
          rasa.r[,itau] <- rasa.beta * rasa.r[,itau] +
            (1 - rasa.beta) * colSums(asl(grad_Theta,itau)^2) / p
          rasa.lhat[,itau] <- pmax(rasa.lhat[,itau], rasa.l[,itau])
          rasa.rhat[,itau] <- pmax(rasa.rhat[,itau], rasa.r[,itau])
        }
      }
      if (it == numIterStage1) step_size <- armijo.alpha
      if (it > numIterStage1) {
        # Stage 2: Adapt smoothness + RASA updates ----------------
        for (itau in 1:ntau) {
          sg_obj <- sg.rasa(
            asl(Theta,itau), sigma2[itau], lambda[,itau],
            asl(grad_Theta,itau), grad_zeta[itau], grad_eta[,itau],
            step_size[itau], step_decay,
            rasa.l[,itau], rasa.r[,itau], rasa.lhat[,itau], rasa.rhat[,itau],
            rasa.beta, consts = list(G = G))
          Theta[,,itau] <- sg_obj$Theta
          sigma2[itau] <- sg_obj$sigma2
          lambda[,itau] <- sg_obj$lambda
          rasa.l[,itau] <- sg_obj$rasa.l; rasa.r[,itau] <- sg_obj$rasa.r
          rasa.lhat[,itau] <- sg_obj$rasa.lhat; rasa.rhat[,itau] <- sg_obj$rasa.rhat
        }
      }
      
      ## log
      # TODO: modify for the case where accumulate_data=TRUE
      if (!accumulate_data) {
        MSE <- N / (N + n) * MSE + n / (N + n) * mse_curr
      } else {
        MSE <- mse_curr
      }
      N <- N + n
      # m_mean <- m_mean * N / (N + n) + m / (N + n)
      # m2_mean <- m2_mean * N / (N + n) + m2 / (N + n)
      
      if (record.estimate) {
        est_trace$Theta <- cbind(est_trace$Theta, c(Theta[,,1]))
        est_trace$theta_mu <- cbind(est_trace$theta_mu, theta_mu)
        est_trace$sigma2 <- c(est_trace$sigma2, sigma2[1])
        est_trace$lambda <- cbind(est_trace$lambda, lambda[,1])
        est_trace$avg_mse <- c(est_trace$avg_mse, MSE)
      }
      
      time.end <- Sys.time()
      if (record.runtime)
        runtime <- c(runtime, as.numeric(time.end - time.start))
      niter <- niter + 1
    }
    
    ## update tuning parameters -------------
    ord <- order(MSE, decreasing = FALSE)
    mse_min_idx <- ord[1:ntau_parent]
    parent_tau_idx <- (mse_min_idx - 1) %/% nexplore + 1
    tau_path <- c(tau_path, list(
      list(tau = tau[mse_min_idx],
           parent = tau_parent[parent_tau_idx],
           parent_id = parent_tau_idx)))
    tau_parent <- tau[mse_min_idx]
    Theta_parent[,,] <- Theta[,,mse_min_idx]
    sigma2_parent[] <- sigma2[mse_min_idx]
    lambda_parent[,] <- lambda[,mse_min_idx]
    step_size <- step_size[mse_min_idx]
    step_size_enlarge_factor <- step_size_enlarge_factor[mse_min_idx]
    stepsize_path <- c(stepsize_path, step_size)
    rasa.l[,] <- rasa.l[,mse_min_idx]
    rasa.r[,] <- rasa.r[,mse_min_idx]
    rasa.lhat[,] <- rasa.lhat[,mse_min_idx]
    rasa.rhat[,] <- rasa.rhat[,mse_min_idx]
    MSE <- MSE[mse_min_idx]
  }
  
  out <- list(
    Theta = Theta, theta_mu = theta_mu, sigma2 = sigma2, lambda = lambda,
    tau_mu = tau_mu, tau = tau, tau_path = tau_path,
    MSE = MSE, stepsize_path = stepsize_path
  )
  
  if (record.estimate) {
    est_trace$Theta <- structure(est_trace$Theta, dim = c(p,q,ncol(est_trace$Theta)))
    out[['est_trace']] <- est_trace
  }
  
  if (record.runtime)
    out[['runtime']] <- runtime
  
  # output
  return(out)
}


fpca.pspc.batch <- function(
    Ly, Lt,
    basis = NULL, init = NULL,
    nbasis = 10, npc = NULL,
    tau_mu = 10^(-4), tau = 10^(0),
    sampling = FALSE, nsample = 100,
    step_size = 1e-3, momentum = 0.9,
    step_size_ratio = "sqinv", step_size_decay_factor = 0.1,
    maxIter = 100, eps.conv = 1e-4,
    verbose = FALSE,
    record.estimate = FALSE,
    record.objvalues = FALSE) {
  
  if (is.null(basis)) {
    stop("NOT IMPLEMENTED")
    init <- NULL
  }
  stopifnot(class(basis) %in% c("basisfd", "TensorBasis"))
  if (is(basis, "basisfd")) {
    p <- basis$nbasis
    is.basisfd <- TRUE
  } else if (is(basis, "TensorBasis")) {
    p <- attr(basis, "nbasis")
    is.basisfd <- FALSE
  }
  
  G <- get_basis_inprod_matrix(basis)
  Omega <- get_basis_penalty_matrix(basis, penLfd = 2)
  
  if (is.null(init) || !("theta_mu" %in% names(init))) {
    # TODO: avoid contructing basis matrix on all points
    smth.mean <- smooth_basis.gcv(
      argvals = bind_argvals(Lt, is.basisfd),
      Y = do.call(c, Ly), basis = basis, lambda = tau_mu)
    theta_mu <- smth.mean$fd$coef
  } else {theta_mu <- init[['theta_mu']]}
  
  if (is.null(init) || !all((c("Theta", "sigma2", "lambda") %in% names(init)))) {
    # TODO: initialization with a subset of data
    stop("NOT IMPLEMENTED")
  }
  Theta <- init[['Theta']]
  sigma2 <- init[['sigma2']]
  lambda <- init[['lambda']]
  
  q <- npc
  if (length(lambda) < q) {
    lambda <- c(lambda, rep(min(lambda) * 0.1, q - length(lambda)))
  } else lambda <- lambda[1:q]
  
  if (is.null(step_size_ratio))
    step_size_ratio <- "sqinv"
  if (is.numeric(step_size_ratio)) {
    stopifnot(length(step_size_ratio) >= q)
    step_size_ratio <- step_size_ratio[1:q]
  } else if (is.character(step_size_ratio)) {
    if (step_size_ratio == "sqinv") {
      step_size_ratio <- lambda^(-2)
    } else if (step_size_ratio == "identity") {
      step_size_ratio <- rep(1, q)
    } else {
      stop("Invalid step_size_ratio")}
  }
  step_size_ratio <- step_size_ratio / step_size_ratio[1]
  
  if (length(tau) > 1) {
    stop("NOT IMPLEMENTED for multiple tau's")
  }
  
  stopifnot(is.list(Lt) && is.list(Ly))
  N <- length(Lt)
  stopifnot(length(Ly) == N)
  
  cumgrad_Theta <- matrix(0, nrow = p, ncol = q)
  cumgrad_zeta <- 0
  cumgrad_eta <- numeric(q)
  MSE <- c()
  objvalues <- c()
  
  if (record.estimate) {
    est_trace <- list(Theta = c(Theta), sigma2 = sigma2,
                      lambda = lambda)
  }
  
  for (it in 1:maxIter) {
    if (verbose) message("Iter ", it)
    
    Theta_old <- Theta 
    sigma2_old <- sigma2
    lambda_old <- lambda
    
    # gradients
    if (!sampling) {
      optObj <- pspc_stats(
        Ly, Lt, basis, Theta, theta_mu, sigma2, lambda,
        stats = c("grad", "mse"))
    } else {
      sub_idx <- ceiling(runif(nsample) * N)
      optObj <- pspc_stats(
        Ly[sub_idx], Lt[sub_idx], basis, Theta, theta_mu, sigma2, lambda,
        stats = c("grad", "mse"))
    }
    
    grad_Theta <- optObj$grad_Theta + tau * pspc_smoothpen_grad(Theta, Omega)
    grad_Theta <- as.matrix(manifold.Stiefel.project(grad_Theta, Theta, G))
    grad_zeta <- optObj$grad_zeta
    grad_eta <- optObj$grad_eta
    mse_curr <- optObj$mse
    MSE <- c(MSE, mse_curr)
    
    if (record.objvalues) {
      if (!sampling) {
        objvalues[it] <- mse_curr + tau * pspc_smoothpen(Theta, Omega)
      } else {
        mse_full <- pspc_stats(
          Ly, Lt, basis, Theta, theta_mu, sigma2, lambda,
          stats = c("mse"))$mse
        objvalues[it] <- mse_full + tau * pspc_smoothpen(Theta, Omega)
      }
    }
    
    if (it == 1) {
      cumgrad_Theta <- grad_Theta
      cumgrad_zeta <- grad_zeta
      cumgrad_eta <- grad_eta
    } else {
      # FIXME: cumgrad has been implemented in sg.vanilla
      # here are repeated calculations
      cumgrad_Theta <- momentum * cumgrad_Theta + (1 - momentum) * grad_Theta
      cumgrad_zeta <- momentum * cumgrad_zeta + (1 - momentum) * grad_zeta
      cumgrad_eta <- momentum * cumgrad_eta + (1 - momentum) * grad_eta
    }
    
    # update parameters
    step_decay <- (1 / it)^step_size_decay_factor
    # sg_obj <- sg.vanilla(
    #   Theta, sigma2, lambda,
    #   grad_Theta, grad_zeta, grad_eta,
    #   cumgrad_Theta, cumgrad_zeta, cumgrad_eta,
    #   step_size, step_decay, momentum,
    #   consts = list(G = G)) 
    # Theta <- sg_obj$Theta 
    # sigma2 <- sg_obj$sigma2
    # lambda <- sg_obj$lambda
    Theta <- Theta - step_size * step_decay *
      sweep(cumgrad_Theta, 2, step_size_ratio, "*")
    sigma2 <- sigma2 * exp(-step_size * step_decay * cumgrad_zeta)
    lambda <- lambda * exp(
      -step_size * step_decay * step_size_ratio * cumgrad_eta)
    Theta <- orthonormalize.basiscoef(Theta, G)
    
    if (record.estimate) {
      est_trace$Theta <- cbind(est_trace$Theta, c(Theta))
      est_trace$sigma2 <- c(est_trace$sigma2, sigma2)
      est_trace$lambda <- cbind(est_trace$lambda, lambda)
    }
    
    if (max(abs(c(Theta - Theta_old, sigma2 - sigma2_old,
                  lambda - lambda_old))) < eps.conv) break
  }
  
  out <- list(
    Theta = Theta, theta_mu = theta_mu, sigma2 = sigma2, lambda = lambda,
    tau_mu = smth.mean$lambda.min, tau = tau, MSE = MSE,
    objvalues = objvalues
  )
  
  if (record.estimate) {
    est_trace$Theta <- array(est_trace$Theta, dim = c(p, q, ncol(est_trace$Theta)))
    out[['est_trace']] <- est_trace
  }
  
  return(out)
}


fpca.pspc.init <- function() {
  
}



