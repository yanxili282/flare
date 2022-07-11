#### Author: Yanxi Li, yli282@g.uky.edu
### Required packages
library(MASS)
library(emg)
library(mixtools)
library(ggplot2)

### load functions:
flare.reg.boot <- function(dat, B, n){
  flare.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    flare.boot[[i]] <-  try.flare.reg(y =dat.boot$y,
                                      x = dat.boot$x, emg = 0)
  }
  return(flare.boot)
}

emg.reg.boot <- function(dat, B, n){
  emg.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    emg.boot[[i]] <-  emg.reg(x = cbind(1, dat.boot$x), 
                              y = dat.boot$y) 
  }
  return(emg.boot)
}

ols.reg.boot <- function(dat, B, n){
  ols.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    ols.boot[[i]] <-  lm(dat.boot$y~dat.boot$x)
  }
  return(ols.boot)
}

regmix.reg.boot <- function(dat, B, n){
  regmix.boot <- list()
  for (i in 1:n){
    dat.boot <- dat[sample(1:nrow(dat), B, replace=TRUE),]
    regmix.boot[[i]] <-  regmixEM(y = dat.boot$y,
                                  x = dat.boot$x)
  }
  return(regmix.boot)
}

## Produce Web Figure 1 in the Supporting Information:
set.seed(12)
n <- 200
x <- rnorm(n, 0, 1)
beta <- c(-2, 4)
sigma <- 0.5
alpha <- 0.05

X <- cbind(1,x)
E <- remg(n, mu = 0, sigma = sigma, lambda = alpha)
Y <- X %*% beta + E
X <- X[,-1]
simdat <- data.frame(X, Y)
ggplot(simdat, aes(x=X, y=Y)) + geom_point(size = 0.8, color = "#00A087B2")+
  xlab("x") + ylab("y")+ 
  geom_abline(intercept = -2, slope = 4, color = "red")



set.seed(10)
### Produce Table 1 and Table 2 in the main paper
n <- 200
x <- rnorm(n, 0, 1)

## Table 1:
beta <- c(-2, 4)
sigma <- 0.5 
alpha <- 0.05

E <- remg(n, mu = 0, sigma = sigma, lambda = alpha)
y <- cbind(1,x) %*% beta + E

dat1 <- data.frame(x, y)

out1 <- emg.reg(x = cbind(1,x), y = y)
out1[3:6]

boot1 <- emg.reg.boot(dat=dat1, B=1000, n=100)

emg_BETA1 <- c()
emg_BETA2 <- c()
emg_SIGMA <- c()
emg_ALPHA <- c()

for (i in 1:100){
  emg_BETA1[i] <- boot1[[i]]$beta[1]
  emg_BETA2[i] <- boot1[[i]]$beta[2]
  emg_SIGMA[i] <- boot1[[i]]$sigma
  emg_ALPHA[i] <- boot1[[i]]$ex.rate
}

sd(emg_BETA1)
sd(emg_BETA2)
sd(emg_SIGMA)
sd(emg_ALPHA)

## Table 2:
set.seed(12)
n <- 200
x <- rnorm(n, 0, 1)

lambda <- c(0.6, 0.4)
beta <- c(-2, 4)
sigma <- 0.5
alpha <- 0.05

E <- c()
Z <- c()
for (i in 1:n){
  indicator <- runif(1, min = 0, max = 1)
  if (indicator <= lambda[1]){
    eps <- rnorm(1, mean=0, sd=sigma)
    z <- "Gaussian"
    E <- c(E, eps)
    Z <- c(Z, z)
  } else{
    eps <- rexp(1, rate=alpha)
    z <- "exponential"
    E <- c(E, eps)
    Z <- c(Z, z)
  }
}

E <- matrix(data = E, nrow = n, ncol = 1, byrow = FALSE,
            dimnames = NULL)

y <- cbind(1, x) %*% beta + E
dat2 <- data.frame(x, y)
out2 <- try.flare.reg(y, x, emg = 0)
out2[4:7]
## standard error by bootstrap
boot2.1 <- flare.reg.boot(dat=dat2, B=1000, n=100)

flare_LAMBDA <- c()
flare_BETA1 <- c()
flare_BETA2 <- c()
flare_SIGMA <- c()
flare_ALPHA <- c()

for (i in 1:100){
  flare_LAMBDA[i] <-boot2.1[[i]]$lambda[1]
  flare_BETA1[i] <- boot2.1[[i]]$beta[1]
  flare_BETA2[i] <- boot2.1[[i]]$beta[2]
  flare_SIGMA[i] <- sqrt(boot2.1[[i]]$sigma)
  flare_ALPHA[i] <- boot2.1[[i]]$alpha
}

sd(flare_LAMBDA[flare_SIGMA<2])
sd(flare_BETA1[flare_SIGMA<2])
sd(flare_BETA2[flare_SIGMA<2])
sd(flare_SIGMA[flare_SIGMA<2])
sd(flare_ALPHA[flare_SIGMA<2])

## standard error by Fisher Information
out2.sd <- sd.flare.reg(Y=y, X=x, lambda = out2$lambda, beta = c(out2$beta),
                        sigma = sqrt(out2$sigma), alpha = out2$alpha)

### Produce Web Figure 2 in the Supporting Information
dat2$Group <- Z

## Count observations from each component
length(dat2$Group[dat2$Group=="Gaussian"])
length(dat2$Group[dat2$Group=="exponential"])

ggplot(dat2, aes(x=x, y=y, colour = Group)) + geom_point(size = 0.8, aes(shape=Group))+
  xlab("x") + ylab("y")

dat2[out2$posterior[,2]>=0.8,]

Cluster <- rep("Gaussian", 200)
Cluster[out2$posterior[,2]>=0.8] <- "exponential"
dat2$Cluster <- Cluster

ggplot(dat2, aes(x=x, y=y, colour = Cluster)) + geom_point(size = 0.8, aes(shape= Cluster))+
  xlab("x") + ylab("y")

## Count observations classified to each component:
length(dat2$Cluster[dat2$Cluster=="Gaussian"])
length(dat2$Cluster[dat2$Cluster=="exponential"])

## Produce Web Table 0 in the Supporting Information
dat2
