linesearch.pck <- function(
    k, Theta, lambda, sigma2, theta_mu, tau,
    objfun, grad_thetak, grad_etak, f0, fprev = NA,
    ls.alpha = 1e-4, ls.beta = 0.5, ls.optimism = 2, ls.maxiter = 25,
    opt.minStepsize = 1e-10, opt.maxStepsize = 1000,
    opt.initStepsize = 1, opt.finalStepsize = 1, opt.accuracy = 1e-6,
    opt.minGradnorm = .Machine$double.eps) {
  # global constant: Ly, Lt, B, Cov, G, Omega
  
  ngf <- sqrt(sum(grad_thetak^2 + grad_etak^2))
  
  if (is.na(fprev) || abs(fprev - f0) < opt.accuracy) {
    stepsize <- opt.initStepsize / ngf
  } else {
    stepsize <- 2 * (fprev - f0) / ngf^2
    # Look a little further
    stepsize <- stepsize * ls.optimism
  }
  
  stepcount <- 1
  ThetaNew <- Theta
  lambdaNew <- lambda
  while (stepcount <= ls.maxiter) {
    ThetaNew[,k] <- Theta[,k] - stepsize * grad_thetak
    if (k == 1) {
      ThetaNew[,k] <- manifold.Stiefel.retractSinglePC(ThetaNew[,k], NULL, G)
    } else {
      ThetaNew[,k] <- manifold.Stiefel.retractSinglePC(
        ThetaNew[,k], ThetaNew[,1:(k-1),drop=F], G)
    }
    lambdaNew[k] <- lambda[k] * exp(-stepsize * grad_etak)
    fnew <- objfun(ThetaNew, lambdaNew, sigma2, theta_mu, tau)$fval
    
    if (fnew > f0 - ls.alpha * stepsize * ngf^2) {
      stepsize <- stepsize * ls.beta
      stepcount <- stepcount + 1
    } else {
      Theta[,k] <- ThetaNew[,k]
      lambda[k] <- lambdaNew[k]
      break
    }
  }
  
  if (stepcount == ls.maxiter) {
    stepsize <- 0
  }
  
  # convergence check
  stepsize.rescale <- stepsize * ngf
  if (stepsize.rescale < opt.minStepsize || ngf < opt.minGradnorm) {
    opt.conv <- TRUE
  } else {
    opt.conv <- FALSE
  }
  
  return(list(
    thetak = Theta[,k], lambdak = lambda[k], opt.conv = opt.conv,
    stepsize = stepsize, stepsize.rescale = stepsize.rescale
  ))
}


linesearch.sigma2 <- function(
    Theta, lambda, sigma2, theta_mu, tau,
    objfun, grad_zeta, f0, fprev = NA,
    ls.alpha = 1e-4, ls.beta = 0.5, ls.optimism = 2, ls.maxiter = 25,
    opt.minStepsize = 1e-10, opt.maxStepsize = 1000,
    opt.initStepsize = 1, opt.finalStepsize = 1, opt.accuracy = 1e-6,
    opt.minGradnorm = .Machine$double.eps) {
  # global constant: Ly, Lt, B, Cov, G, Omega
  
  ngf <- abs(grad_zeta)
  
  if (is.na(fprev) || abs(fprev - f0) < opt.accuracy) {
    stepsize <- opt.initStepsize / ngf
  } else {
    stepsize <- 2 * (fprev - f0) / ngf^2
    # Look a little further
    stepsize <- stepsize * ls.optimism
  }
  
  stepcount <- 1
  while (stepcount <= ls.maxiter) {
    sigma2New <- sigma2 * exp(-stepsize * grad_zeta)
    fnew <- objfun(Theta, lambda, sigma2New, theta_mu, tau)$fval
    
    if (fnew > f0 - ls.alpha * stepsize * ngf^2) {
      stepsize <- stepsize * ls.beta
      stepcount <- stepcount + 1
    } else {
      sigma2 <- sigma2New
      break
    }
  }
  
  if (stepcount == ls.maxiter) {
    stepsize <- 0
  }
  
  # convergence check
  stepsize.rescale <- stepsize * ngf
  if (stepsize.rescale < opt.minStepsize || ngf < opt.minGradnorm) {
    opt.conv <- TRUE
  } else {
    opt.conv <- FALSE
  }
  
  return(list(
    sigma2 = sigma2, opt.conv = opt.conv,
    stepsize = stepsize, stepsize.rescale = stepsize.rescale
  ))
}