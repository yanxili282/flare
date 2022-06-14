#### Author: Yanxi Li, yli282@g.uky.edu
### Required packages:
library(MASS)
library(emg)
library(mixtools)

## load functions:
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
  # print(start)
  mle(function(mu, sigma, lambda) {
    emg.nllik(x, mu, sigma, lambda)
  }, method = "L-BFGS-B", lower = lower, upper = upper, start = start, 
  ...)
}
######################################################################
emg.reg <- function(x, y, beta=NULL, sigma=NULL, lambda=NULL, maxit = 500, epsilon = 1e-04){
  ## y is a vector, x is a matrix with intercept
  ## here our lambda is alpha in the flare setting
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
######################################################################
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
      if (any(class(beta.new) == "try-error") || sum(is.na(beta.new)) > 
          0){
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
      #print(c(beta.new, lambda.new, sigma.new, alpha.new))
      diff <- c(lambda.new, beta.new, sigma.new, alpha.new) - 
        c(lambda, beta, sigma, alpha)
      z.new2 = Z(res.new, sigma.new, lambda.new, alpha.new)
      Q.new <- Q(res.new, sigma.new, lambda.new, alpha.new, 
                 z.new2)
      q.diff = Q.new - Q1
      #print (c(Q.new, Q1))
      if (q.diff < 0) 
        m = m + 1
      else m = 101
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
######################################################################
######################################################################
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
######################################################################
BIC.flare <- function(out){
  n <- length(out$y)
  L <- out$loglik
  k <- length(out$beta)+3
  BIC <- k*log(n)-2*L
  return (BIC)
}
######################################################################
######################################################################
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
######################################################################
BIC.emg <- function(out){
  n <- length(out$residual)
  L <- out$loglik
  k <- length(out$beta)+2
  BIC <- k*log(n)-2*L
  return (BIC)
}
######################################################################
######################################################################
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

######################################################################
######################################################################
BIC.regmixEM <- function(out){
  n <- length(out$y)
  L <- out$loglik
  k <- length(out$beta[,1])*2+3
  BIC <- k*log(n)-2*L
  return (BIC)
}
######################################################################
######################################################################
######################################################################
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

### Perform large simulation for setting M1:
### Produce row 1--3 for Web Table 1:
## n=100
set.seed(201)
B <- 1000

M1_list <- list()
for (i in 1:B){
  M1_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/20)
}
outM1_100 <- run_on_sim(B=B, M_list=M1_list)

outM1_100
write.table(outM1_100[[1]], file="outM1_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(201)
B <- 1000
M1_list <- list()
for (i in 1:B){
  M1_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/20)
}
outM1_500 <- run_on_sim(B=B, M_list=M1_list)

outM1_500

write.table(outM1_500[[1]], file="outM1_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(201)
B <- 1000
M1_list <- list()
for (i in 1:B){
  M1_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/20)
}
outM1_1000 <- run_on_sim(B=B, M_list=M1_list)

outM1_1000

write.table(outM1_1000[[1]], file="outM1_1000.txt", row.names=TRUE, col.names=TRUE)

### Perform large simulation for setting M2:
### Produce row 4--6 for Web Table 1:
## n=100
set.seed(201)
B <- 1000
M2_list <- list()
for (i in 1:B){
  M2_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/6)
}
outM2_100 <- run_on_sim(B=B, M_list=M2_list)

outM2_100

write.table(outM2_100[[1]], file="outM2_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(201)
B <- 1000
M2_list <- list()
for (i in 1:B){
  M2_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/6)
}
outM2_500 <- run_on_sim(B=B, M_list=M2_list)

outM2_500

write.table(outM2_500[[1]], file="outM2_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(201)
B <- 1000
M2_list <- list()
for (i in 1:B){
  M2_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/6)
}
outM2_1000 <- run_on_sim(B=B, M_list=M2_list)

outM2_1000

### Perform large simulation for setting M3:
### Produce row 7--9 for Web Table 1:
## n=100
set.seed(201)
B <- 1000
M3_list <- list()
for (i in 1:B){
  M3_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/2)
}
outM3_100 <- run_on_sim(B=B, M_list=M3_list)

outM3_100

write.table(outM3_100[[1]], file="outM3_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(201)
B <- 1000
M3_list <- list()
for (i in 1:B){
  M3_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/2)
}
outM3_500 <- run_on_sim(B=B, M_list=M3_list)

outM3_500

write.table(outM3_500[[1]], file="outM3_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(201)
B <- 1000
M3_list <- list()
for (i in 1:B){
  M3_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                 lambda=c(1/3, 2/3), sigma=1/2, alpha=1/2)
}
outM3_1000 <- run_on_sim(B=B, M_list=M3_list)

outM3_1000

write.table(outM3_1000[[1]], file="outM3_1000.txt", row.names=TRUE, col.names=TRUE)

### Perform large simulation for setting M4:
### Produce row 10--12 for Web Table 1:
## n=150
set.seed(204)
B <- 1000
M4_list <- list()
for (i in 1:B){
  M4_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=150, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/20)
}
outM4_150 <- run_on_sim(B=B, M_list=M4_list)

outM4_150

write.table(outM4_150[[1]], file="outM4_150.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(204)
B <- 1000
M4_list <- list()
for (i in 1:B){
  M4_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/20)
}
outM4_500 <- run_on_sim(B=B, M_list=M4_list)

outM4_500

write.table(outM4_500[[1]], file="outM4_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(203)
B <- 1000
M4_list <- list()
for (i in 1:B){
  M4_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/20)
}
outM4_1000 <- run_on_sim(B=B, M_list=M4_list)

outM4_1000

write.table(outM4_1000[[1]], file="outM4_1000.txt", row.names=TRUE, col.names=TRUE)

### Perform large simulation for setting M5:
### Produce row 13--15 for Web Table 1:
## n=100
set.seed(23)
B <- 1000
M5_list <- list()
for (i in 1:B){
  M5_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/6)
}
outM5_100 <- run_on_sim(B=B, M_list=M5_list)

outM5_100

write.table(outM5_100[[1]], file="outM5_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(20)
B <- 1000
M5_list <- list()
for (i in 1:B){
  M5_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/6)
}
outM5_500 <- run_on_sim(B=B, M_list=M5_list)

outM5_500

write.table(outM5_500[[1]], file="outM5_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(20)
B <- 1000
M5_list <- list()
for (i in 1:B){
  M5_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/6)
}
outM5_1000 <- run_on_sim(B=B, M_list=M5_list)

outM5_1000

write.table(outM5_1000[[1]], file="outM5_1000.txt", row.names=TRUE, col.names=TRUE)

### Perform large simulation for setting M6:
### Produce row 16--18 for Web Table 1:
## n=100
set.seed(23)
B <- 1000
M6_list <- list()
for (i in 1:B){
  M6_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)
}
outM6_100 <- run_on_sim(B=B, M_list=M6_list)

outM6_100

write.table(outM6_100[[1]], file="outM6_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(23)
B <- 1000
M6_list <- list()
for (i in 1:B){
  M6_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)
}
outM6_500 <- run_on_sim(B=B, M_list=M6_list)

outM6_500

write.table(outM6_500[[1]], file="outM6_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(23)
B <- 1000
M6_list <- list()
for (i in 1:B){
  M6_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                 lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)
}
outM6_1000 <- run_on_sim(B=B, M_list=M6_list)

outM6_1000

write.table(outM6_1000[[1]], file="outM6_1000.txt", row.names=TRUE, col.names=TRUE)


### Perform large simulation for setting M7:
### Produce row 19--21 for Web Table 1:
## n=100
set.seed(201)
B <- 1000
M7_list <- list()
for (i in 1:B){
  M7_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/2, alpha=1/25)
}
outM7_100 <- run_on_sim(B=B, M_list=M7_list)

outM7_100

write.table(outM7_100[[1]], file="outM7_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(201)
B <- 1000
M7_list <- list()
for (i in 1:B){
  M7_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/2, alpha=1/25)
}
outM7_500 <- run_on_sim(B=B, M_list=M7_list)

outM7_500

write.table(outM7_500[[1]], file="outM7_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(201)
B <- 1000
M7_list <- list()
for (i in 1:B){
  M7_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/2, alpha=1/25)
}
outM7_1000 <- run_on_sim(B=B, M_list=M7_list)

outM7_1000

write.table(outM7_1000[[1]], file="outM7_1000.txt", row.names=TRUE, col.names=TRUE)


### Perform large simulation for setting M8:
### Produce row 22--24 for Web Table 1:
## n=100
set.seed(202)
B <- 1000
M8_list <- list()
for (i in 1:B){
  M8_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/3, alpha=1/5)
}
outM8_100 <- run_on_sim(B=B, M_list=M8_list)

outM8_100

write.table(outM8_100[[1]], file="outM8_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(202)
B <- 1000
M8_list <- list()
for (i in 1:B){
  M8_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/3, alpha=1/5)
}
outM8_500 <- run_on_sim(B=B, M_list=M8_list)

outM8_500

write.table(outM8_500[[1]], file="outM8_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(202)
B <- 1000
M8_list <- list()
for (i in 1:B){
  M8_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/3, alpha=1/5)
}
outM8_1000 <- run_on_sim(B=B, M_list=M8_list)

outM8_1000

write.table(outM8_1000[[1]], file="outM8_1000.txt", row.names=TRUE, col.names=TRUE)

### Perform large simulation for setting M9:
### Produce row 25--27 for Web Table 1:
## n=100
set.seed(201)
B <- 1000
M9_list <- list()
for (i in 1:B){
  M9_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/2, alpha=1/2)
}
outM9_100 <- run_on_sim(B=B, M_list=M9_list)

outM9_100

write.table(outM9_100[[1]], file="outM9_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(201)
B <- 1000
M9_list <- list()
for (i in 1:B){
  M9_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/2, alpha=1/2)
}
outM9_500 <- run_on_sim(B=B, M_list=M9_list)

outM9_500

write.table(outM9_500[[1]], file="outM9_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(201)
B <- 1000
M9_list <- list()
for (i in 1:B){
  M9_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                 lambda=c(1/2, 1/2), sigma=1/2, alpha=1/2)
}
outM9_1000 <- run_on_sim(B=B, M_list=M9_list)

outM9_1000

write.table(outM9_1000[[1]], file="outM9_1000.txt", row.names=TRUE, col.names=TRUE)

### Perform large simulation for setting M10:
### Produce row 28--30 for Web Table 1:
## n=100
set.seed(23)
B <- 1000
M10_list <- list()
for (i in 1:B){
  M10_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/25)
}
outM10_100 <- run_on_sim(B=B, M_list=M10_list)

outM10_100

write.table(outM10_100[[1]], file="outM10_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(23)
B <- 1000
M10_list <- list()
for (i in 1:B){
  M10_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/25)
}
outM10_500 <- run_on_sim(B=B, M_list=M10_list)

outM10_500

write.table(outM10_500[[1]], file="outM10_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(23)
B <- 1000
M10_list <- list()
for (i in 1:B){
  M10_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/25)
}
outM10_1000 <- run_on_sim(B=B, M_list=M10_list)

outM10_1000

write.table(outM10_1000[[1]], file="outM10_1000.txt", row.names=TRUE, col.names=TRUE)


### Perform large simulation for setting M11:
### Produce row 31--33 for Web Table 1:
## n=100
set.seed(25)
B <- 1000
M11_list <- list()
for (i in 1:B){
  M11_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=0.2)
}
outM11_100 <- run_on_sim(B=B, M_list=M11_list)

outM11_100

write.table(outM11_100[[1]], file="outM11_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(25)
B <- 1000
M11_list <- list()
for (i in 1:B){
  M11_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=0.2)
}
outM11_500 <- run_on_sim(B=B, M_list=M11_list)

outM11_500

write.table(outM11_500[[1]], file="outM11_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(23)
B <- 1000
M11_list <- list()
for (i in 1:B){
  M11_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=0.2)
}
outM11_1000 <- run_on_sim(B=B, M_list=M11_list)

outM11_1000

write.table(outM11_1000[[1]], file="outM11_1000.txt", row.names=TRUE, col.names=TRUE)


### Perform large simulation for setting M12:
### Produce row 34--36 for Web Table 1:
## n=100
set.seed(23)
B <- 1000
M12_list <- list()
for (i in 1:B){
  M12_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=0.5, alpha=0.5)
}
outM12_100 <- run_on_sim(B=B, M_list=M12_list)

outM12_100

write.table(outM12_100[[1]], file="outM12_100.txt", row.names=TRUE, col.names=TRUE)

## n=500
set.seed(23)
B <- 1000
M12_list <- list()
for (i in 1:B){
  M12_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=0.5, alpha=0.5)
}
outM12_500 <- run_on_sim(B=B, M_list=M12_list)

outM12_500

write.table(outM12_500[[1]], file="outM12_500.txt", row.names=TRUE, col.names=TRUE)

## n=1000
set.seed(23)
B <- 1000
M12_list <- list()
for (i in 1:B){
  M12_list[[i]] <- sim.flare_reg (xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=0.5, alpha=0.5)
}
outM12_1000 <- run_on_sim(B=B, M_list=M12_list)

outM12_1000

write.table(outM12_1000[[1]], file="outM12_1000.txt", row.names=TRUE, col.names=TRUE)

