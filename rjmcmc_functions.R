# Libraries
library(MASS)
library(bspec)

# log likelihood
llike <- function(yt, X, a, sigma2.eps) {
  k = length(a)  # Same as ncol(X)
  n = length(yt)
  if (k == 1) {
    ll = -n * log(sqrt(2 * pi * sigma2.eps)) -0.5 * sum((yt - X * a) ^ 2) / sigma2.eps
  }
  else {
    ll = -n * log(sqrt(2 * pi * sigma2.eps)) -0.5 * sum((yt - X %*% a) ^ 2) / sigma2.eps
  }
  return(ll)
}

# Joint log prior
lprior <- function(a, sigma2.eps, sigma2.a, alpha.0, beta.0) {
  k = length(a)
  lprior.a = -k * log(sqrt(2 * pi * sigma2.a)) - sum(a ^ 2) / (2 * sigma2.a) 
  lprior.sigma2.eps = -(alpha.0 + 1) * log(sigma2.eps) - beta.0 / sigma2.eps
  lp = lprior.a + lprior.sigma2.eps
  return(lp)
}

# Joint log posterior
lpost <- function(yt, X, a, sigma2.eps, sigma2.a, alpha.0, beta.0) {
  lpst = lprior(a, sigma2.eps, sigma2.a, alpha.0, beta.0) + 
    llike(yt, X, a, sigma2.eps)
  return(lpst)
}