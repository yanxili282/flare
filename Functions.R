#### Author: Yanxi Li, yli282@g.uky.edu
### Required packages
library(emg)
library(MASS)
library(mixtools)
library(brms)
library(gamlss)

### function to estimate the EMG regression:
## y is a vector, x is a matrix with intercept (first column is vector 1)
emg.reg <- function(x, y, beta=NULL, sigma=NULL, lambda=NULL, maxit = 10000, epsilon = 1e-04){      
  emg.mle2 <- function (x, lower = NULL, upper = NULL, start = NULL, ...) 
  {
    if (is.null(lower)) {
      lower <- list(mu = min(x), sigma = sd(x)/10, lambda = abs(0.001/mean(x)))
    }
    if (is.null(upper)) {
      upper <- list(mu = max(x), sigma = (max(x) - min(x))/4, 
                    lambda = abs(100/mean(x)))
    }
    if (is.null(start)) {
      start <- list(mu = mean(x) - sd(x) * ((skewness(x)/2)^(1/3)), 
                    sigma = sqrt(abs(sd(x)^2 * (1 - (skewness(x)/2)^(2/3)))), 
                    lambda = (1/(sd(x) * abs((skewness(x))/2)^(1/3))))
      if (is.nan(start$mu)) {
        start <- list(mu = (lower$mu + upper$mu)/2, sigma = (lower$sigma + 
                                                               upper$sigma)/2, lambda = abs((lower$lambda + upper$lambda)/2))
      }
    }
    mle(function(mu, sigma, lambda) {
      emg.nllik(x, mu, sigma, lambda)
    }, method = "L-BFGS-B", lower = lower, upper = upper, start = start, 
    ...)
  }
  y <- matrix(y)
  N <- length(y)
  if (is.null(beta)) {
    lm.out = lsfit(x, y ,intercept=FALSE)
    beta = lm.out$coef
    beta = as.numeric(c(beta))
  }
  if (is.null(sigma)) {
    lm.out = lm(y ~ x[, 2])
    sigma = rexp(1, rate = sqrt(1/anova(lm.out)$Mean[2]))
  }
  if (is.null(lambda)) {
    lm.out = lm(y ~ x[, 2])
    a = 1/sum(lm.out$res[lm.out$res > 0])
    lambda = abs(rnorm(1, a))
  }
  mu = 0
  beta_log_fn <- function(beta){
    emg.nllik(c(y-x%*%beta), mu = mu, sigma = sigma, lambda =lambda)
  }
  sum.finite <- function(x) {
    sum(x[is.finite(x)])
  }
  Q <- 0
  Q[2] <- -emg.nllik(c(y-x%*%beta), mu = mu, sigma = sigma, lambda =lambda)
  k <- 2
  iter <- 0
  while (abs(Q[k]-Q[k-1])>=epsilon && iter < maxit){
    iter = iter + 1
    result <- optim(par = beta, beta_log_fn)
    beta <- result$par
    out <- emg.mle2(c(y-x%*%beta))
    mu <- 0
    sigma <- as.numeric(coef(out)[2]) 
    lambda <- as.numeric(coef(out)[3]) 
    beta_log_fn <- function(beta){
      E <- y-x%*%beta
      emg.nllik(c(E), mu = mu, sigma = sigma, lambda = lambda)
    }
    k <- k + 1
    Q[k] <- -emg.nllik(c(y-x%*%beta), mu = mu, sigma = sigma, lambda =lambda)
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  fit <- list("all.loglik" = Q,
              "residual" =c(y-x%*%beta),
              "loglik" = Q[k],
              "beta"=beta,
              "sigma"=sigma,
              "ex.rate"=lambda)
  return (fit)
}

### function to estimate the flare regression:
## y is a vector, x is a matrix without intercept (first column is NOT vector 1)
try.flare.reg <- function (y, x, lambda = NULL, beta = NULL, sigma = NULL, alpha = NULL, emg = 1,
                           epsilon = 1e-04, maxit = 10000, verb = FALSE, restart = 50)
{
  lambda = lambda[1]
  loglik <- function(res, sigma, lambda, alpha) {
    tmp <- lambda * dnorm(res, sd = sqrt(sigma)) + (1 - 
                                                      lambda) * dexp(res, rate = alpha)
    sum(log(tmp))
  }
  Q = function(res, sigma, lambda, alpha, z) {
    Q <- sum(z * log(lambda)) + sum((1 - z) * log(1 - lambda)) - 
      log(2 * pi * sigma) * sum(z)/2 - sum(z * res^2)/2/sigma + 
      log(alpha) * sum(1 - z) - alpha * sum((1 - z) * 
                                              res)
    Q
  }
  Z <- function(res, sigma, lambda, alpha) {
    z = rep(1, length(res))
    z[res > 0] = lambda/(lambda + (1 - lambda) * sqrt(2 * 
                                                        pi * sigma) * alpha * as.vector(exp(res[res > 0]^2/2/sigma - 
                                                                                              alpha * res[res > 0])))
    z
  }
  x <- cbind(1, x)
  n <- length(y)
  p <- ncol(x)
  if (is.null(beta)) {
    lm.out = lsfit(x, y ,intercept=FALSE)
    beta = lm.out$coef
    beta = as.numeric(c(beta))
  }
  if (is.null(sigma)) {
    lm.out = lm(y ~ x[, 2])
    sigma = rexp(1, rate = sqrt(1/anova(lm.out)$Mean[2]))
  }
  if (is.null(alpha)) {
    lm.out = lm(y ~ x[, 2])
    a = 1/sum(lm.out$res[lm.out$res > 0])
    alpha = abs(rnorm(1, a))
  }
  if (emg == 1){
    est = try(emg.reg(x, y, beta = beta, sigma = sigma, lambda = alpha), silent = TRUE)
    beta = est$beta
    sigma = est$sigma
    alpha = est$ex.rate
  } else{
    beta = beta
    sigma = sigma
    alpha = alpha
  }
  if (is.null(lambda)){
    mem <- kmeans(y,2)$cluster
    y1 <- y[mem==1]
    y2 <- y[mem==2]
    lambda1 <- length(y1)/length(y)
    lambda2 <- 1-lambda1
    lambda <- lambda1
  }
  diff <- 1
  iter <- 0
  counts <- 0
  ll.counts <- 0
  xbeta <- x %*% beta
  res <- y - xbeta
  dn <- dnorm(res, sd = sqrt(sigma))
  de <- dexp(res, rate = alpha)
  obsloglik <- loglik(res, sigma, lambda, alpha)
  ll <- obsloglik
  Q1 <- -Inf
  all.Q <- NULL
  z = Z(res, sigma, lambda, alpha)
  while (sum(abs(diff) > epsilon) > 0 && iter < maxit) {
    iter = iter + 1
    temp = (ginv(-1/sigma * t(x) %*% sweep(x, 1, z, "*")) %*% (1/sigma * (apply(sweep(x, 1, z * 
                                                        (y - x %*% beta), "*"), 2, sum)) + alpha * apply(sweep(x, 
                                                              1, 1 - z, "*"), 2, sum)))
    m = 1
    while (m < restart){
      beta.new <- try(beta - temp, silent = TRUE)
      if (any(class(beta.new) == "try-error") || sum(is.na(beta.new)) >0){
        Q_beta <- function(beta_hat){
          res_hat <- y - x%*% beta_hat
          Q(res_hat, sigma, lambda, alpha, z)
        }
        result <- optim(par = beta, Q_beta)
        beta.new <- result$par
      } else{
        beta.new = beta.new
      }
      xbeta.new <- x %*% beta.new
      res.new <- y - xbeta.new
      Q.beta <- Q(res.new, sigma, lambda, alpha, z)
      z.new = Z(res.new, sigma, lambda, alpha)
      lambda.new <- mean(z.new)
      sigma.new <- sum(z.new * (res.new^2))/sum(z.new)
      alpha.new <- sum(1 - z.new[res.new > 0])/sum((1 - 
                                                      z.new[res.new > 0]) * res.new[res.new > 0])
      diff <- c(lambda.new, beta.new, sigma.new, alpha.new) - 
        c(lambda, beta, sigma, alpha)
      z.new2 = Z(res.new, sigma.new, lambda.new, alpha.new)
      Q.new <- Q(res.new, sigma.new, lambda.new, alpha.new, 
                 z.new2)
      q.diff = Q.new - Q1
      if (q.diff < 0) 
        m = m + 1
      else m = 101
    }
    if (m == restart){
      print("Too many attempts at step-halving!")
    }
    lambda <- lambda.new
    beta <- beta.new
    xbeta <- xbeta.new
    res <- res.new
    sigma <- sigma.new
    alpha <- alpha.new
    z <- z.new2
    newobsloglik <- loglik(res.new, sigma.new, lambda.new, 
                           alpha.new)
    ll <- c(ll, newobsloglik)
    counts <- counts + (Q.new < Q1)
    all.Q <- c(all.Q, Q.new)
    ll.counts <- ll.counts + (obsloglik > newobsloglik)
    Q1 <- Q.new
    obsloglik <- newobsloglik
    if (verb == TRUE) 
      cat("iteration=", iter, "diff=", diff, "log-likelihood", 
          obsloglik, "\n")
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  a = list(x = x, y = y, posterior = cbind(z, 1 - z), lambda = c(lambda, 
                                                                 1 - lambda), beta = beta, sigma = sigma, alpha = alpha, 
           loglik = obsloglik, all.loglik = ll, ft = "flaremixEM")
  class(a) = "mixEM"
  a
}  

### function to estimate standard errors for point estimates from try.flare.reg()(Louis's method)
## Y is a vector, X is a matrix without intercept (first column is NOT vector 1)
## lambda, beta, sigma, and alpha are point estimates for their corresponding estimates
## lambda is a vector, others are scalars
sd.flare.reg <- function(Y, X, lambda, beta, sigma, alpha){
  X <- cbind(1, X)
  m <- length(beta)
  m2 <- m +3
  N <- length(Y)
  Z <- Y - X%*%beta
  gamma_est <- function(z, lambda, sigma, alpha){
    gamma_sum <- lambda[1]*dnorm(z, mean = 0, sd = sigma, log = FALSE) + lambda[2]*dexp(z, rate = alpha, log = FALSE)
    if (gamma_sum >0){
      gamma1 <- lambda[1]*dnorm(z, mean = 0, sd = sigma, log = FALSE)/gamma_sum
      gamma2 <- lambda[2]*dexp(z, rate = alpha, log = FALSE)/gamma_sum
    } else{
      gamma1 <- 1
      gamma2 <- 0
    }
    return (c(gamma1, gamma2))
  }
  partial_l_lambda_lambda <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <- -gamma_est(z, lambda, sigma, alpha)[1]/(lambda[1]^2) - gamma_est(z, lambda, sigma, alpha)[2]/(lambda[2]^2)
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  partial_l_sigma_sigma <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]*(1/(2*sigma^4)-z^2/(sigma^6))
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  partial_l_alpha_alpha <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <- -gamma_est(z, lambda, sigma, alpha)[2]/(alpha^2)
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  partial_l_beta_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- -gamma_est(z, lambda, sigma, alpha)[1]*x^2/(sigma^2)
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  partial_l_sigma_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- -gamma_est(z, lambda, sigma, alpha)[1]*(x*z/(sigma^4))
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  partial_l_alpha_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[2]*x
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  partial_l_betaj_betak <- function(j, k, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x1 <- X[i,j]
      x2 <- X[i,k]
      z <- c(Z[i,])
      partial <- -gamma_est(z, lambda, sigma, alpha)[1]*x1*x2/(sigma^2)
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  I1 <- matrix(0, nrow = m2, ncol = m2)
  I1[1, 1] <- partial_l_lambda_lambda(Z, lambda, sigma, alpha)
  I1[2, 2] <- partial_l_sigma_sigma(Z, lambda, sigma, alpha)
  I1[3, 3] <- partial_l_alpha_alpha(Z, lambda, sigma, alpha)
  for (j in 4:m2){
    I1[j, j] <- partial_l_beta_beta(j-3, Y, X, beta, lambda, sigma, alpha)
  }
  for (j in 2:m2){
    I1[1, j] <- 0
  }
  I1[2, 3] <- 0
  for (k in 4:m2){
    I1[2, k] <- partial_l_sigma_beta(k-3, Y, X, beta, lambda, sigma, alpha)
  }
  for (k in 4:m2){
    I1[3, k] <- partial_l_alpha_beta(k-3, Y, X, beta, lambda, sigma, alpha)
  }
  for (j in 4:(m2-1)){
    for (k in (j+1):m2){
      I1[j, k] <- partial_l_betaj_betak(j-3, k-3, Y, X, beta, lambda, sigma, alpha)
    }
  }
  I1 <- I1 + t(I1)
  E_l_lambda_l_lambda <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]/(lambda[1]^2) +gamma_est(z, lambda, sigma, alpha)[2]/(lambda[2]^2)
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_lambda_l_sigma <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <-  gamma_est(z, lambda, sigma, alpha)[1]*((-1/(2*sigma^2))+z^2/(2*sigma^4))/lambda[1]
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_lambda_l_alpha <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <--gamma_est(z, lambda, sigma, alpha)[2]*((1/alpha)-z)/lambda[2]
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_lambda_l_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]*(z*x/sigma^2)/lambda[1]-gamma_est(z, lambda, sigma, alpha)[2]*alpha*x/lambda[2]
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_sigma_l_sigma <- function(Z, lambda, sigma, alpha){
    partial_sum <- 0
    N <- length(Z)
    for (i in 1:N){
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]*(-1/(2*sigma^2)+z^2/(2*sigma^4))^2
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_sigma_l_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]*(-1/(2*sigma^2)+z^2/(2*sigma^4))*(z*x/(sigma^2))
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_alpha_l_alpha <- function(Z, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[2]*(1/alpha - z)^2
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_alpha_l_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[2]*(z/alpha)*alpha*x
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_beta_beta <- function(j, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x <- X[i,j]
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]*(z*x/(sigma^2)) + gamma_est(z, lambda, sigma, alpha)[2]*(alpha*x)^2
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  E_l_betaj_betak <- function(j, k, Y, X, beta, lambda, sigma, alpha){
    Z <- Y - X%*%beta
    N <- length(Z)
    partial_sum <- 0
    for (i in 1:N){
      x1 <- X[i,j]
      x2 <- X[i,k]
      z <- c(Z[i,])
      partial <- gamma_est(z, lambda, sigma, alpha)[1]*z^2*x1*x2/(sigma^4)+gamma_est(z, lambda, sigma, alpha)[2]*alpha^2*x1*x2
      partial_sum <- partial_sum + partial
    }
    return (partial_sum)
  }
  I2 <- matrix(0, nrow = m2, ncol = m2)
  I2[1, 1] <- E_l_lambda_l_lambda(Z, lambda, sigma, alpha)
  I2[2, 2] <- E_l_sigma_l_sigma(Z, lambda, sigma, alpha)
  I2[3, 3] <- E_l_alpha_l_alpha(Z, lambda, sigma, alpha)
  for (j in 4:m2){
    I2[j, j] <- E_l_beta_beta(j-3, Y, X, beta, lambda, sigma, alpha)
  }
  I2[1, 2] <- E_l_lambda_l_sigma(Z, lambda, sigma, alpha)
  I2[1, 3] <- E_l_lambda_l_alpha(Z, lambda, sigma, alpha)
  for (j in 4:m2){
    I2[1, j] <- E_l_lambda_l_beta(j-3, Y, X, beta, lambda, sigma, alpha)
  }
  I2[2, 3] <- 0
  for (j in 4:m2){
    I2[2, j] <- E_l_sigma_l_beta(j-3, Y, X, beta, lambda, sigma, alpha)
  }
  for (j in 4:m2){
    I2[3, j] <- E_l_alpha_l_beta(j-3, Y, X, beta, lambda, sigma, alpha)
  }
  for (j in 4:(m2-1)){
    for (k in (j+1):m2){
      I2[j, k] <- E_l_betaj_betak(j-3, k-3, Y, X, beta, lambda, sigma, alpha)
    }
  }
  I2 <- I2 + t(I2)
  partial_li_theta <- function(i, m, Y, X, beta, lambda, sigma, alpha){
    partial_mat <- matrix(NA, nrow = m+3, ncol = 1)
    Z <- Y - X%*%beta
    z <- c(Z[i,])
    partial_mat[1, 1] <- gamma_est(z, lambda, sigma, alpha)[1]/lambda[1] - gamma_est(z, lambda, sigma, alpha)[2]/lambda[2]
    partial_mat[2, 1] <- gamma_est(z, lambda, sigma, alpha)[1]*(-1/(2*sigma^2)+z^2/(2*sigma^4))
    partial_mat[3, 1] <- gamma_est(z, lambda, sigma, alpha)[2]*(1/alpha - z)
    for (j in 1:m){
      x <- X[i,j]
      partial_mat[j+3, 1] <- gamma_est(z, lambda, sigma, alpha)[1]*(z*x/(sigma^2))+gamma_est(z, lambda, sigma, alpha)[2]*x*alpha
    }
    return (partial_mat)
  }
  mat <- partial_li_theta(N, m, Y, X, beta, lambda, sigma, alpha)%*%(t(partial_li_theta(i=1, m, Y, X, beta, lambda, sigma, alpha)))
  I3 <- matrix(0, nrow = m2, ncol = m2)
  matrix_sum_i <- function(i, m, Y, X, beta, lambda, sigma, alpha){
    m2 <- m+3
    mat_sum <- matrix(0, nrow = m2, ncol = m2)
    for (k in (1+i):N){
      mat <- partial_li_theta(i, m, Y, X, beta, lambda, sigma, alpha)%*%(t(partial_li_theta(k, m, Y, X, beta, lambda, sigma, alpha)))
      mat_sum <- mat_sum + mat
    }
    return (mat_sum)
  }
  for (i in 1:(N-1)){
    I3 <- I3 + matrix_sum_i(i, m, Y, X, beta, lambda, sigma, alpha)
  }
  Z <- Y - X%*%beta
  mat_sum <-  matrix(0, nrow = m2, ncol = 1)
  for (i in 1:N){
    z <- c(Z[i,])
    mat <- matrix(NA, nrow=m2, ncol=1)
    mat[1, 1] <- gamma_est(z, lambda, sigma, alpha)[1]/lambda[1] - gamma_est(z, lambda, sigma, alpha)[2]/lambda[2]
    mat[2, 1] <- gamma_est(z, lambda, sigma, alpha)[1]*(-1/(2*sigma^2)+z^2/(2*sigma^4))
    mat[3, 1] <- gamma_est(z, lambda, sigma, alpha)[2]*(1/alpha - z)
    for (j in 1:m){
      x <- X[i,j]
      mat[j+3, 1] <- gamma_est(z, lambda, sigma, alpha)[1]*(z*x/(sigma^2)) +gamma_est(z, lambda, sigma, alpha)[2]*x*alpha
    }
    mat_sum <- mat_sum + mat
  }
  I4 <- mat_sum%*%(t(mat_sum)) 
  I <- -I1-I2-2*I3 + I4
  est_sd <- sqrt(diag(ginv(I)))
  output <- list("The estimated covariance matrix" = ginv(I),
                 "The estimated standard deviation for lambda" =est_sd[1],
                 "The estimated standard deviation for sigma^2"=est_sd[2],
                 "The estimated standard deviation for alpha"=est_sd[3],
                 "The estimated standard deviation for beta"=est_sd[4:m2])
  return(output)
}

### function to estimate standard errors for point estimates from try.flare.reg()(Bootstrap method)
## dat is the input dataframe; B is the bootstrapping sample size
flare.reg.boot <- function(dat, B, n){
  flare.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    flare.boot[[i]] <-  try.flare.reg(y =dat.boot$y,
                                      x = dat.boot$x, emg = 0)
  }
  return(flare.boot)
}

### function to estimate standard errors for point estimates from emg.reg()(Bootstrap method)
## dat is the input dataframe; B is the bootstrapping sample size
emg.reg.boot <- function(dat, B, n){
  emg.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    emg.boot[[i]] <-  emg.reg(x = cbind(1, dat.boot$x), 
                              y = dat.boot$y) 
  }
  return(emg.boot)
}

### function to estimate standard errors for point estimates from ls()(Bootstrap method)
## dat is the input dataframe; B is the bootstrapping sample size
ols.reg.boot <- function(dat, B, n){
  ols.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    ols.boot[[i]] <-  lm(dat.boot$y~dat.boot$x)
  }
  return(ols.boot)
}

### function to estimate standard errors for point estimates from regmix()(Bootstrap method)
## dat is the input dataframe; B is the bootstrapping sample size
regmix.reg.boot <- function(dat, B, n){
  regmix.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    regmix.boot[[i]] <-  regmixEM(y = dat.boot$y,
                                  x = dat.boot$x)
  }
  return(regmix.boot)
}

### function to simulate n dataframes assuming the flare regression model
## the output covariate matrix x is a matrix without intercept
sim.flare_reg <- function(xmin, xmax, n, beta, lambda, sigma, alpha){
  m = length(beta)
  X <- matrix(runif(n*(m-1), min = xmin, max = xmax), 
              nrow = n, ncol = (m-1), byrow = FALSE, dimnames = NULL)
  X <- cbind(1,X)
  lambda1 = lambda[1]
  
  E <- c()
  for (i in 1:n){
    indicator <- runif(1, min = 0, max = 1)
    if (indicator <= lambda1){
      eps <- rnorm(1, mean=0, sd=sigma)
      E <- c(E, eps)
    } else{
      eps <- rexp(1, rate=alpha)
      E <- c(E, eps)
    }
  }
  
  E <- matrix(data = E, nrow = n, ncol = 1, byrow = FALSE,
              dimnames = NULL)
  
  Y <- X %*% beta + E
  X <- X[,-1]
  sim.out <- list("X" = X,
                  "Y" = Y)
  return(sim.out)
}

### function to simulate n dataframes assuming the EMG regression model
## the output covariate matrix x is a matrix without intercept
sim.emg_reg <- function(xmin, xmax, n, beta, sigma, alpha){
  m = length(beta)
  X <- matrix(runif(n*(m-1), min = xmin, max = xmax), 
              nrow = n, ncol = (m-1), byrow = FALSE, dimnames = NULL)
  X <- cbind(1,X)
  E <- remg(n, mu = 0, sigma = sigma, lambda = alpha)
  Y <- X %*% beta + E
  X <- X[,-1]
  sim.out <- list("X" = X,
                  "Y" = Y)
  return(sim.out)
}

### function to simulate n dataframes assuming the mixture of linear regressions model with two components
## the output covariate matrix x is a matrix without intercept
sim.regmix <- function(xmin, xmax, n, lambda, beta1, beta2, sigma1, sigma2){
  m = length(beta1)
  X <- matrix(runif(n*(m-1), min = xmin, max = xmax), 
              nrow = n, ncol = (m-1), byrow = FALSE, dimnames = NULL)
  X <- cbind(1,X)
  
  lambda1 = lambda[1]
  
  Y <- c()
  for (i in 1:n){
    indicator <- runif(1, min = 0, max = 1)
    if (indicator <= lambda1){
      eps <- rnorm(1, mean=0, sd=sigma1)
      y <- X[i,] %*% beta1 + eps
      Y <- c(Y, c(y))
    } else{
      eps <- rnorm(1, mean=0, sd=sigma2)
      y <- X[i,] %*% beta2 + eps
      Y <- c(Y, c(y))
    }
  }
  X <- X[,-1]
  sim.out <- list("X" = X,
                  "Y" = Y)
  return(sim.out)
}

### function to calculate BIC value for the fitted flare regression model
## out is the output from function try.flare()
BIC.flare <- function(out){
  n <- length(out$y)
  L <- out$loglik
  k <- length(out$beta)+3
  BIC <- k*log(n)-2*L
  return (BIC)
}

### function to calculate BIC value for the fitted EMG regression model
## out is the output from function emg.reg()
BIC.emg <- function(out){
  n <- length(out$residual)
  L <- out$loglik
  k <- length(out$beta)+2
  BIC <- k*log(n)-2*L
  return (BIC)
}

### function to calculate BIC value for the fitted mixture of linear regressions model with two components
## out is the output from function regmixEM()
BIC.regmixEM <- function(out){
  n <- length(out$y)
  L <- out$loglik
  k <- length(out$beta[,1])*2+3
  BIC <- k*log(n)-2*L
  return (BIC)
}

### modified version for the function regmixEM() in the package "mixtools"
regmixEM2 <- function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2, 
                       addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
                       maxit = 10000, verb = FALSE) 
{
  if (arbmean == FALSE && arbvar == FALSE) {
    stop(paste("Must change constraints on beta and/or sigma!", 
               "\n"))
  }
  s = sigma
  if (addintercept) {
    x = cbind(1, x)
  }
  n <- length(y)
  p <- ncol(x)
  tmp <- regmix.init(y = y, x = x, lambda = lambda, beta = beta, 
                     s = s, k = k, addintercept = addintercept, arbmean = arbmean, 
                     arbvar = arbvar)
  lambda <- tmp$lambda
  beta <- tmp$beta
  s <- tmp$s
  k <- tmp$k
  diff <- 1
  iter <- 0
  xbeta <- x %*% beta
  res <- (y - xbeta)^2
  if (arbmean == FALSE) {
    res <- sapply(1:k, function(i) res)
  }
  comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
                                                             s^2)))))
  obsloglik <- sum(log(apply(comp, 1, sum)))
  ll <- obsloglik
  z = matrix(nrow = n, ncol = k)
  restarts <- 0
  while (diff > epsilon && iter < maxit) {
    for (i in 1:n) {
      for (j in 1:k) {
        z.denom = c()
        for (h in 1:k) {
          z.denom = c(z.denom, (lambda[h]/lambda[j]) * 
                        (s[j * arbvar + (1 - arbvar)]/s[h * arbvar + 
                                                          (1 - arbvar)]) * exp(-0.5 * ((1/s[h * 
                                                                                              arbvar + (1 - arbvar)]^2) * res[i, h] - 
                                                                                         (1/s[j * arbvar + (1 - arbvar)]^2) * res[i, 
                                                                                                                                  j])))
        }
        z[i, j] = 1/sum(z.denom)
      }
    }
    z = z/apply(z, 1, sum)
    lambda.new <- apply(z, 2, mean)
    if (sum(lambda.new < 1e-08) > 0 || is.na(sum(lambda.new))) {
      sing <- 1
    }
    else {
      if (arbmean == FALSE) {
        if (addintercept) {
          beta.new <- lm(y ~ x[, -1], weights = apply(t(t(z)/(s^2)), 
                                                      1, sum))$coef
        }
        else beta.new <- lm(y ~ x - 1, weights = apply(t(t(z)/(s^2)), 
                                                       1, sum))$coef
      }
      else {
        if (addintercept) {
          lm.out <- lapply(1:k, function(i) lm(y ~ x[, 
                                                     -1], weights = z[, i]))
        }
        else lm.out <- lapply(1:k, function(i) lm(y ~ 
                                                    x - 1, weights = z[, i]))
        beta.new <- sapply(lm.out, coef)
      }
      xbeta.new <- x %*% beta.new
      res <- (y - xbeta.new)^2
      if (arbmean == FALSE) {
        res <- sapply(1:k, function(i) res)
      }
      if (arbvar) {
        s.new <- sqrt(sapply(1:k, function(i) sum(z[, 
                                                    i] * (res[, i]))/sum(z[, i])))
      }
      else s.new <- sqrt(sum(z * res)/n)
      lambda <- lambda.new
      beta <- beta.new
      xbeta <- x %*% beta
      s <- s.new
      sing <- sum(s < 1e-08)
      comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, 
                                                        xbeta[, i * arbmean + (1 - arbmean)], s[i * 
                                                                                                  arbvar + (1 - arbvar)]))
      comp <- sapply(comp, cbind)
      compsum <- apply(comp, 1, sum)
      newobsloglik <- sum(log(compsum))
    }
    if (sing > 0 || is.na(newobsloglik) || newobsloglik < 
        obsloglik || abs(newobsloglik) == Inf) {
      restarts <- restarts + 1
      if (restarts > 15) 
        tmp <- regmix.init(y = y, x = x, k = k, addintercept = addintercept, 
                           arbmean = arbmean, arbvar = arbvar)
      lambda <- tmp$lambda
      beta <- tmp$beta
      s <- tmp$s
      k <- tmp$k
      diff <- 1
      iter <- 0
      xbeta <- x %*% beta
      res <- (y - xbeta)^2
      if (arbmean == FALSE) {
        res <- sapply(1:k, function(i) res)
      }
      comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
                                                                 s^2)))))
      obsloglik <- sum(log(apply(comp, 1, sum)))
      ll <- obsloglik
    }
    else {
      diff <- newobsloglik - obsloglik
      obsloglik <- newobsloglik
      ll <- c(ll, obsloglik)
      iter <- iter + 1
      if (verb) {
        cat("iteration=", iter, "diff=", diff, "log-likelihood", 
            obsloglik, "\n")
      }
    }
  }
  scale.order = order(s)
  sigma.min = min(s)
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  if (arbmean == FALSE) {
    z = z[, scale.order]
    names(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a = list(x = x, y = y, lambda = lambda[scale.order], 
             beta = beta, sigma = sigma.min, scale = s[scale.order]/sigma.min, 
             loglik = obsloglik, posterior = z[, scale.order], 
             all.loglik = ll, restarts = restarts, ft = "regmixEM")
    class(a) = "mixEM"
    a
  }
  else {
    rownames(beta) <- c(paste("beta", ".", 0:(p - 1), sep = ""))
    colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a = list(x = x, y = y, lambda = lambda, beta = beta, 
             sigma = s, loglik = obsloglik, posterior = z, all.loglik = ll, 
             restarts = restarts, ft = "regmixEM")
    class(a) = "mixEM"
    a
  }
}


### function to perform large simulation study:
## assuming B is the bootstrap number, M is simulated data
## assuming our simulated X is without the intercept
run_on_sim <- function(B, M_list){
  M.flare_list <- list()
  M.emg_list <- list()
  M.lm_list <- list()
  M.regmix_list <- list()
  b <- 1
  while(b <=B){
    M.flare_list[[b]] <- try(try.flare.reg(y = M_list[[b]]$Y, x = M_list[[b]]$X, emg=0))
    if ("try-error" %in% class(M.flare_list[[b]])){
      M.flare_list[[b]] <- try(try.flare.reg(y = M_list[[b]]$Y, x = M_list[[b]]$X, emg=1))
    }else if(M.flare_list[[b]]$sigma>=10){
      M.flare_list[[b]] <- try(try.flare.reg(y = M_list[[b]]$Y, x = M_list[[b]]$X, emg=1))
    }else{
      M.flare_list[[b]] = M.flare_list[[b]] 
      M.flare_list[[b]] = append(M.flare_list[[b]][4:7], list("BIC"=BIC.flare(M.flare_list[[b]])))
      res <-  c(M_list[[b]]$Y-(cbind(1, M_list[[b]]$X)%*%M.flare_list[[b]]$beta+M.flare_list[[b]]$lambda[2]/(M.flare_list[[b]]$alpha)))
      mse <- sum(res^2)/length(res)
      M.flare_list[[b]] = append(M.flare_list[[b]], list("MSE"=mse))
      b = b+1
    }
  }
  b <- 1
  while(b <= B){
    M.emg_list[[b]] <- try(emg.reg(x= cbind(1, M_list[[b]]$X), y=M_list[[b]]$Y))
    if ("try-error" %in% class(M.emg_list[[b]])){
      M.emg_list[[b]] <- try(emg.reg(x= cbind(1, M_list[[b]]$X), y=M_list[[b]]$Y), beta=rnorm(length(c(M.flare_list[[1]]$beta))))
    } else{
      M.emg_list[[b]] = M.emg_list[[b]]
      M.emg_list[[b]] = append(M.emg_list[[b]][3:6], list("BIC"=BIC.emg(M.emg_list[[b]])))
      res <-  c(M_list[[b]]$Y-cbind(1, M_list[[b]]$X)%*%M.emg_list[[b]]$beta-1/M.emg_list[[b]]$ex.rate)
      mse <- sum(res^2)/length(res)
      M.emg_list[[b]] = append(M.emg_list[[b]], list("MSE"=mse))
      b = b+1
    }
  }
  b <- 1
  while(b <= B){
    M.lm_list[[b]] = lm(M_list[[b]]$Y~M_list[[b]]$X)
    res <- c(M.lm_list[[b]]$resid)
    mse <- sum(res^2)/length(res)
    M.lm_list[[b]] = append((as.numeric(M.lm_list[[b]]$coefficients)), list("sigma"=summary(M.lm_list[[b]])$sigma))
    M.lm_list[[b]] =append(M.lm_list[[b]], list("BIC"=BIC(lm(M_list[[b]]$Y~M_list[[b]]$X))))
    M.lm_list[[b]] = append(M.lm_list[[b]], list("MSE"=mse))
    b = b+1
  }
  b <- 1
  while(b <=B){
    M.regmix_list[[b]] = regmixEM2(y=c(M_list[[b]]$Y), x=M_list[[b]]$X, addintercept = TRUE)
    M.regmix_list[[b]] =append(M.regmix_list[[b]][3:6], list("BIC"=BIC.regmixEM(M.regmix_list[[b]])))
    res <- c(M_list[[b]]$Y - M.regmix_list[[b]]$lambda[1]* cbind(1, M_list[[b]]$X)%*%M.regmix_list[[b]]$beta[,1]
             - M.regmix_list[[b]]$lambda[2]* cbind(1, M_list[[b]]$X)%*%M.regmix_list[[b]]$beta[,2])
    mse <- sum(res^2)/length(res)
    M.regmix_list[[b]] =append(M.regmix_list[[b]], list("MSE"=mse))
    b = b+1
  }
  minBIC_index <- c()
  for(i in 1:B){
    BIC_vec <- c(M.flare_list[[i]]$BIC, M.emg_list[[i]]$BIC, M.lm_list[[i]]$BIC, M.regmix_list[[i]]$BIC)
    minBIC <- min(BIC_vec)
    minBIC_index[i] <- match(minBIC, BIC_vec)
  }
  percent1 <- sum(minBIC_index == 1)/B
  percent2 <- sum(minBIC_index == 2)/B
  percent3 <- sum(minBIC_index == 3)/B
  percent4 <- sum(minBIC_index == 4)/B
  C <- length(c(unlist(M.flare_list[[1]]), unlist(M.emg_list[[1]]), 
                unlist(M.lm_list[[1]]), unlist(M.regmix_list[[1]])))
  C1 <- length(unlist(M.flare_list[[1]]))
  C2 <- length(unlist(M.emg_list[[1]]))
  C3 <- length(unlist(M.lm_list[[1]]))
  C4 <- length(unlist(M.regmix_list[[1]]))
  sim_vec <- c()
  for (j in 1:B){
    temp_vec <- c(unlist(M.flare_list[[j]]), unlist(M.emg_list[[j]]), 
                  unlist(M.lm_list[[j]]), unlist(M.regmix_list[[j]]))
    sim_vec <- c(sim_vec, temp_vec)
  }
  sim_mat <- matrix(sim_vec, nrow = B, ncol = C, byrow = TRUE)
  beta_length <- length(c(M.flare_list[[1]]$beta))
  C1_name <- rep("flare", C1)
  C2_name <- rep("exGaussian", C2)
  C3_name <- rep("reg", C3)
  C4_name <- rep("regmix", C4)
  C1_name[1:2] <- c("prop1 (flare)", "prop2 (flare)")
  C1_name[3:(3+beta_length-1)] <- rep("beta (flare)", beta_length)
  C1_name[3+beta_length] <- "sigma (flare)"
  C1_name[3+beta_length+1] <- "ex rate (flare)"
  C1_name[3+beta_length+2] <- "BIC (flare)"
  C1_name[3+beta_length+3] <- "MSE (flare)"
  C2_name[1] <- "loglik (emg)"
  C2_name[2:(2+beta_length-1)] <- rep("beta (emg)", beta_length)
  C2_name[2+beta_length] <- "sigma (emg)"
  C2_name[2+beta_length+1] <- "ex rate (emg)"
  C2_name[2+beta_length+2] <- "BIC (emg)"
  C2_name[2+beta_length+3] <- "MSE (emg)"
  C3_name[1:beta_length] <- rep("beta (reg)", beta_length)
  C3_name[beta_length+1] <- "sigma (reg)"
  C3_name[beta_length+2] <- "BIC (reg)"
  C3_name[beta_length+3] <- "MSE (reg)"
  C4_name[1:2] <- c("prop1 (regmix)", "prop2 (regmix)")
  C4_name[3:(3+beta_length-1)] <- rep("beta1 (regmix)", beta_length)
  C4_name[(3+beta_length):(3+2*beta_length-1)] <- rep("beta2 (regmix)", beta_length)
  C4_name[(3+2*beta_length):(3+2*beta_length+1)] <- c("sigma1 (regmix)", "sigma2 (regmix)")
  C4_name[3+2*beta_length+2] <- "loglik (regmix)"
  C4_name[3+2*beta_length+3] <- "BIC (regmix)"
  C4_name[3+2*beta_length+4] <- "MSE (regmix)"
  col_name <- c(C1_name, C2_name, C3_name, C4_name)
  colnames(sim_mat) <- col_name 
  output <- list("output matrix"=sim_mat,
                 "proportion of flare with lowest BIC"=percent1,
                 "proportion of exGaussian with lowest BIC"=percent2,
                 "proportion of single regression with lowest BIC"=percent3,
                 "proportion of mixed regression with lowest BIC"=percent4)
  output
}

### function to calculate RMSEs and biases from the large simulation study:
## mat is the output from function run_on_sim()
## lambda, beta, sigma, alpha are true parameters
rMSE_bias <- function(mat, lambda, beta, sigma, alpha){
  k <- length(beta)
  flare.lambda <-as.numeric(mat[,1])
  rmse_flare.lambda <- sqrt(mean((flare.lambda-lambda[1])^2))
  bias_flare.lambda <- mean((flare.lambda-lambda[1]))
  rmse_flare.beta <- c()
  bias_flare.beta <- c()
  for (i in 1:k){
    flare.beta <- as.numeric(mat[,(2+i)])
    rmse_flare.beta[i] <- sqrt(mean((flare.beta-beta[i])^2))
    bias_flare.beta[i] <- mean((flare.beta-beta[i]))
  }
  flare.sigma <- sqrt(as.numeric(mat[,"sigma..flare."]))
  rmse_flare.sigma <- sqrt(mean((flare.sigma-sigma)^2))
  bias_flare.sigma <- mean((flare.sigma-sigma))
  flare.alpha <- as.numeric(mat[,"ex.rate..flare."])
  rmse_flare.alpha <- sqrt(mean((flare.alpha-alpha)^2))
  bias_flare.alpha <- mean((flare.alpha-alpha))
  output <- list("rmse"=c(rmse_flare.lambda, rmse_flare.beta, 
                          rmse_flare.sigma, rmse_flare.alpha),
                 "bias"=c(bias_flare.lambda, bias_flare.beta, 
                          bias_flare.sigma, bias_flare.alpha))
  output
}

