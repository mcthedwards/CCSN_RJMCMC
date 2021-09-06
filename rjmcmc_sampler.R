rjmcmc = function(yt, 
                  PC,
                  Niter,
                  burnin,
                  thin = 1,
                  sigma2.a = 1e3,  # Tune this
                  alpha.0 = 1e-3,
                  beta.0 = 1e-3,
                  printerval = 100) {
  
  n = length(yt)  # Sample size
  
  pcN = ncol(PC)  # Number of PCs
  index = 1:pcN   # PC Index
  
  # Standardise PCs
  for (i in 1:pcN) {
    PC[, i] = PC[, i] / sqrt(var(PC[, i]))
  }
  
  # Coerce to matrix
  X = as.matrix(PC)
  
  # Open empty vectors/matrices
  Mp = sigma2.eps = rep(NA, Niter)                # Number of parameters in model
  a = matrix(0, ncol = pcN, nrow = Niter)         # Regression coefficients
  PC.bool = matrix(NA, nrow = Niter, ncol = pcN)  # TRUE/FALSE for which PCs in model
  
  # Initial values
  Mp[1] = pcN                                     # Start with most complex model
  PC.bool[1, ] = index %in% sample(index, Mp[1])  # Populate with TRUE/FALSE
  a[1, PC.bool[1, ]] = rnorm(Mp[1])               # Random starting values
  sigma2.eps[1] = 1                               # Initial variance. Make random?
  
  # Time taken
  ptime = proc.time()[1]
  
  # RJMCMC Sampler
  for (i in 2:Niter) {
    
    #####
    # Print interval
    #####
    if (i %% printerval == 0) {
      print(paste("Iteration", i, ",", "Time elapsed",
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes", ", Mp = ", Mp[i - 1]))
    }

    # Decide what move to make: jump up, jump down or stay
    r <- runif(1)
    if (Mp[i - 1] == 1) {  # Cannot jump down
      if (r < 1 / 2) {
        Mp.star = Mp[i - 1] + 1  # Jump up
      }
      else {
        Mp.star = Mp[i - 1]  # Stay
      }
    }
    else if (Mp[i - 1] == pcN) {  # Cannot jump up
      if (r < 1 / 2) {
        Mp.star = Mp[i - 1] - 1  # Jump down
      }
      else {
        Mp.star = Mp[i - 1]  # Stay
      }
    }
    else {
      if (r < 1 / 3) {
        Mp.star = Mp[i - 1] - 1  # Jump down
      }
      else if (r > 2 / 3)  {
        Mp.star = Mp[i - 1] + 1  # Jump up
      }
      else {
        Mp.star = Mp[i - 1]  # Stay
      }
    }
    
    # Decide the rij's i.e., probability that a proposed jump to model j is attempted from theta_i
    if (Mp[i - 1] == 1) {
      if (Mp.star == 1) {
        rij = rji = 0.5
      }
      if (Mp.star == 2) {
        rij = 0.5
        rji = 1 / 3
      }
    }
    else if (Mp[i - 1] == pcN) {
      if (Mp.star == pcN) {
        rij = rji = 0.5
      }
      if (Mp.star == (pcN - 1)) {
        rij = 0.5
        rji = 1 / 3
      }
    }
    else {
      if (Mp.star == 1) {
        rij = 1 / 3
        rji = 0.5
      }
      else if (Mp.star == pcN) {
        rij = 1 / 3
        rji = 0.5
      }
      else {
        rij = rji = 1 / 3
      }
    }
    
    # Stay: Conjugate priors mean we can sample directly
    if (Mp.star == Mp[i - 1]) {  
      
      PC.star = PC.bool[i - 1, ]  # Old boolean for PC index
      XMp = X[, PC.star]          # Explanatory variables for model
      
      if (Mp.star == 1) {
        XMp = as.matrix(XMp)  # Coerce to matrix if Mp = 1
      }
      
      # Compute the posterior mean and covariance parameters
      XX = t(XMp) %*% XMp
      gXX = ginv(XX)
      beta.hat = gXX %*% t(XMp) %*% yt
      lambda.n = XX + diag(1 / sigma2.a, Mp[i - 1])  # Posterior precision hyperparameter
      glambda.n = ginv(lambda.n)                     # Posterior covariance hyperparameter
      mu.n = glambda.n %*% XX %*% beta.hat           # Posterior mean hyperparameter
      alpha.n = alpha.0 + n / 2
      beta.n = beta.0 + 0.5 * (sum(yt ^ 2) - t(mu.n) %*% lambda.n %*% mu.n)
      
      # Sample parameters
      PC.bool[i, ] = PC.bool[i - 1, ]
      a[i, PC.star] =  c(mvrnorm(1, mu = mu.n, Sigma = glambda.n))    
      sigma2.eps[i] = 1 / rgamma(1, alpha.n, beta.n)
      Mp[i] = Mp[i - 1]
      
    }
    
    # Jump up: Use zeroth order method for automatic proposal scaling
    if (Mp.star > Mp[i - 1]) {  
      
      sigma2 = sigma2.a * (rij / rji) ^ 2  # Using zeroth order method
      
      u = rnorm(1)           # Standard normal random variable
      v = sqrt(sigma2) * u   # Proposal for new coefficient
      # Note we are using "zeroth order method" from Brooks et al. (2003)
      
      # Previous PCs and coefficients
      PC.old = PC.bool[i - 1, ]  # Old boolean for PC index
      a.old = a[i - 1, PC.old]   # Old parameter vector
      
      # Choose PC index to birth
      birth = sample(index[!(index %in% index[PC.old])], 1)
      
      # Proposed PCs (with birth) and coefficients
      PC.star = PC.old          # Call old boolean for PC index
      PC.star[birth] = TRUE     # Attach TRUE to correct index
      a.star = a[i - 1, ]       # Call all previous coefficients
      a.star[birth] = v         # Attach proposed v to correct index
      a.star = a.star[PC.star]  # Keep only the proposed PCs
      
      # Create proposed design matrix
      X.old = X[, PC.old]
      X.star = X[, PC.star]
      
      # log acceptance ratio
      log.alpha = lpost(yt, X.star, a.star, sigma2.eps[i - 1], sigma2.a, alpha.0, beta.0) -
        lpost(yt, X.old, a.old, sigma2.eps[i - 1], sigma2.a, alpha.0, beta.0) +
        log(rji) - log(rij) + log(sqrt(sigma2)) + u ^ 2 / 2
      
      if (log(runif(1)) < log.alpha) {  # ACCEPT
        Mp[i] = Mp.star
        a[i, PC.star] = a.star  # Populate non-zero coefficients in correct index
        sigma2.eps[i] = sigma2.eps[i - 1]
        PC.bool[i, ] = PC.star
      }
      else {  # REJECT
        Mp[i] = Mp[i - 1]
        a[i, ] = a[i - 1, ]
        sigma2.eps[i] = sigma2.eps[i - 1]
        PC.bool[i, ] = PC.old
      }
    }
    
    # Jump down: remove end coefficient
    if (Mp.star < Mp[i - 1]) { 
      
      sigma2 = sigma2.a * (rij / rji) ^ 2   # Zeroth order proposal and jacobian
      
      # Previous PCs and coefficients
      PC.old = PC.bool[i - 1, ]  # Old boolean for PC index
      a.old = a[i - 1, PC.old]   # Old parameter vector
      
      # Choose PC index to kill from current
      death = sample(index[PC.old], 1)
      
      # Find v and u from death
      v = a[i - 1, death]   # Coefficient we drop
      u = v / sqrt(sigma2)  # Rescale 
      
      # Proposed PCs (with birth) and coefficients
      PC.star = PC.old          # Call old boolean for PC index
      PC.star[death] = FALSE    # Attach FALSE to correct index
      a.star = a[i - 1, ]       # Call all previous coefficients
      a.star[death] = 0         # Attach 0 to proposed death index
      a.star = a.star[PC.star]  # Keep only proposed PCs
      
      # Create proposed design matrix
      X.old = X[, PC.old]
      X.star = X[, PC.star]
      
      # log acceptance ratio:  Only negate the proposal/jacobian ratio
      log.alpha = lpost(yt, X.star, a.star, sigma2.eps[i - 1], sigma2.a, alpha.0, beta.0) -
        lpost(yt, X.old, a.old, sigma2.eps[i - 1], sigma2.a, alpha.0, beta.0) +
        log(rji) - log(rij) - log(sqrt(sigma2)) - u ^ 2 / 2
      
      if (log(runif(1)) < log.alpha) {  # ACCEPT
        Mp[i] = Mp.star
        a[i, PC.star] = a.star
        sigma2.eps[i] = sigma2.eps[i - 1]
        PC.bool[i, ] = PC.star
      }
      else {  # REJECT
        Mp[i] = Mp[i - 1]
        a[i, ] = a[i - 1, ]
        sigma2.eps[i] = sigma2.eps[i - 1]
        PC.bool[i, ] = PC.old
      } 
    }
  }
  
  # Post-processing
  keep = seq(burnin, Niter, thin)
  a = a[keep, ]
  sigma2.eps = sigma2.eps[keep]
  Mp = Mp[keep]
  PC.bool = PC.bool[keep, ]
  
  # Reconstruct the signal
  recon <- matrix(NA, nrow = nrow(a), ncol = n)
  for (i in 1:nrow(a)) {
    recon[i, ] = X %*% a[i, ]
  }
  
  # Summarise posterior of signal
  recon.mean = apply(recon, 2, mean)
  recon.05 = apply(recon, 2, quantile, probs = 0.05)
  recon.95 = apply(recon, 2, quantile, probs = 0.95)
  
  # Compute proportion of time each PC is in the MCMC samples
  PC.prop = apply(PC.bool, 2, sum) / length(keep)
  
  # Output 
  output = list(recon.mean = recon.mean,
                recon.05 = recon.05,
                recon.95 = recon.95,
                a = a,
                sigma2.eps = sigma2.eps,
                Mp = Mp,
                PC.bool = PC.bool,
                PC.prop = PC.prop)
  
  class(output) = "rjmcmc"  # Assign S3 class to object
  
  return(output)  # Return output
  
}