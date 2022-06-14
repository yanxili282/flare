#### Author: Yanxi Li, yli282@g.uky.edu
### This manuscript is based on outputs produced from Simulation3.R

## load functions:
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

### Produce Web Table 2 (first 9 rows) and Web Table 3 (first 9 rows)
### Produce Web Figures 7--12 (a) and Web Figures 15--18 (a)
matout <-  read.table(file = "outM1_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M1_100_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.05)

matout <-  read.table(file = "outM1_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M1_500_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.05)

matout <-  read.table(file = "outM1_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M1_1000_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                            sigma = 0.5, alpha = 0.05)

matout <-  read.table(file = "outM2_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M2_100_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                           sigma = 0.5, alpha = 1/6)


matout <-  read.table(file = "outM2_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M2_500_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                           sigma = 0.5, alpha = 1/6)

matout <-  read.table(file = "outM2_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M2_1000_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                            sigma = 0.5, alpha = 1/6)

matout <-  read.table(file = "outM3_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M3_100_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM3_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M3_500_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM3_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M3_1000_output <- rMSE_bias(matout, lambda = c(1/3, 2/3), beta = c(9, 3), 
                            sigma = 0.5, alpha = 0.5)


M123_flare.lambda <- data.frame(sep=rep(c("M3", "M2", "M1"), each=3),
                                size=rep(c(100, 500, 1000),3), 
                                rmse=c(M3_100_output$rmse[1], M3_500_output$rmse[1], M3_1000_output$rmse[1],
                                       M2_100_output$rmse[1], M2_500_output$rmse[1], M2_1000_output$rmse[1],
                                       M1_100_output$rmse[1], M1_500_output$rmse[1], M1_1000_output$rmse[1]),
                                bias=c(M3_100_output$bias[1], M3_500_output$bias[1], M3_1000_output$bias[1],
                                       M2_100_output$bias[1], M2_500_output$bias[1], M2_1000_output$bias[1],
                                       M1_100_output$bias[1], M1_500_output$bias[1], M1_1000_output$bias[1]))

ggplot(M123_flare.lambda, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of lambda with different sample sizes")

ggplot(M123_flare.lambda, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of lambda with sample size")



M123_flare.beta1 <- data.frame(sep=rep(c("M3", "M2", "M1"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M3_100_output$rmse[2], M3_500_output$rmse[2], M3_1000_output$rmse[2],
                                      M2_100_output$rmse[2], M2_500_output$rmse[2], M2_1000_output$rmse[2],
                                      M1_100_output$rmse[2], M1_500_output$rmse[2], M1_1000_output$rmse[2]),
                               bias=c(M3_100_output$bias[2], M3_500_output$bias[2], M3_1000_output$bias[2],
                                      M2_100_output$bias[2], M2_500_output$bias[2], M2_1000_output$bias[2],
                                      M1_100_output$bias[2], M1_500_output$bias[2], M1_1000_output$bias[2]))

ggplot(M123_flare.beta1, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta0 with different sample sizes")

ggplot(M123_flare.beta1, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta0 with different sample sizes")

M123_flare.beta2 <- data.frame(sep=rep(c("M3", "M2", "M1"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M3_100_output$rmse[3], M3_500_output$rmse[3], M3_1000_output$rmse[3],
                                      M2_100_output$rmse[3], M2_500_output$rmse[3], M2_1000_output$rmse[3],
                                      M1_100_output$rmse[3], M1_500_output$rmse[3], M1_1000_output$rmse[3]),
                               bias=c(M3_100_output$bias[3], M3_500_output$bias[3], M3_1000_output$bias[3],
                                      M2_100_output$bias[3], M2_500_output$bias[3], M2_1000_output$bias[3],
                                      M1_100_output$bias[3], M1_500_output$bias[3], M1_1000_output$bias[3]))

ggplot(M123_flare.beta2, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta1 with different sample sizes")

ggplot(M123_flare.beta2, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta1 with different sample sizes")

M123_flare.sigma <- data.frame(sep=rep(c("M3", "M2", "M1"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M3_100_output$rmse[4], M3_500_output$rmse[4], M3_1000_output$rmse[4],
                                      M2_100_output$rmse[4], M2_500_output$rmse[4], M2_1000_output$rmse[4],
                                      M1_100_output$rmse[4], M1_500_output$rmse[4], M1_1000_output$rmse[4]),
                               bias=c(M3_100_output$bias[4], M3_500_output$bias[4], M3_1000_output$bias[4],
                                      M2_100_output$bias[4], M2_500_output$bias[4], M2_1000_output$bias[4],
                                      M1_100_output$bias[4], M1_500_output$bias[4], M1_1000_output$bias[4]))

ggplot(M123_flare.sigma, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of sigma with different sample sizes")

ggplot(M123_flare.sigma, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of sigma with different sample sizes")

M123_flare.alpha <- data.frame(sep=rep(c("M3", "M2", "M1"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M3_100_output$rmse[5], M3_500_output$rmse[5], M3_1000_output$rmse[5],
                                      M2_100_output$rmse[5], M2_500_output$rmse[5], M2_1000_output$rmse[5],
                                      M1_100_output$rmse[5], M1_500_output$rmse[5], M1_1000_output$rmse[5]),
                               bias=c(M3_100_output$bias[5], M3_500_output$bias[5], M3_1000_output$bias[5],
                                      M2_100_output$bias[5], M2_500_output$bias[5], M2_1000_output$bias[5],
                                      M1_100_output$bias[5], M1_500_output$bias[5], M1_1000_output$bias[5]))

ggplot(M123_flare.alpha, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of alpha with different sample sizes")

ggplot(M123_flare.alpha, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of alpha with different sample sizes")


### Produce Web Table 2 (second 9 rows) and Web Table 3 (second 9 rows)
### Produce Web Figures 7--12 (b) and Web Figures 15--18 (b)
matout <-  read.table(file = "outM4_150.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M4_150_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.05)


matout <-  read.table(file = "outM4_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M4_500_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.05)


matout <-  read.table(file = "outM4_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M4_1000_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                            sigma = 0.5, alpha = 0.05)


matout <-  read.table(file = "outM5_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M5_100_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                           sigma = 0.5, alpha = 1/6)


matout <-  read.table(file = "outM5_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M5_500_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                           sigma = 0.5, alpha = 1/6)


matout <-  read.table(file = "outM5_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M5_1000_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                            sigma = 0.5, alpha = 1/6)


matout <-  read.table(file = "outM6_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M6_100_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM6_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M6_500_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                           sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM6_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M6_1000_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(9, 3), 
                            sigma = 0.5, alpha = 0.5)


M456_flare.lambda <- data.frame(sep=rep(c("M6", "M5", "M4"), each=3),
                                size=rep(c(100, 500, 1000),3), 
                                rmse=c(M6_100_output$rmse[1], M6_500_output$rmse[1], M6_1000_output$rmse[1],
                                       M5_100_output$rmse[1], M5_500_output$rmse[1], M5_1000_output$rmse[1],
                                       M4_150_output$rmse[1], M4_500_output$rmse[1], M4_1000_output$rmse[1]),
                                bias=c(M6_100_output$bias[1], M6_500_output$bias[1], M6_1000_output$bias[1],
                                       M5_100_output$bias[1], M5_500_output$bias[1], M5_1000_output$bias[1],
                                       M4_150_output$bias[1], M4_500_output$bias[1], M4_1000_output$bias[1]))

ggplot(M456_flare.lambda, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of lambda with different sample sizes")

ggplot(M456_flare.lambda, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of lambda with different sample sizes")


M456_flare.beta1 <- data.frame(sep=rep(c("M6", "M5", "M4"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M6_100_output$rmse[2], M6_500_output$rmse[2], M6_1000_output$rmse[2],
                                      M5_100_output$rmse[2], M5_500_output$rmse[2], M5_1000_output$rmse[2],
                                      M4_150_output$rmse[2], M4_500_output$rmse[2], M4_1000_output$rmse[2]),
                               bias=c(M6_100_output$bias[2], M6_500_output$bias[2], M6_1000_output$bias[2],
                                      M5_100_output$bias[2], M5_500_output$bias[2], M5_1000_output$bias[2],
                                      M4_150_output$bias[2], M4_500_output$bias[2], M4_1000_output$bias[2]))

ggplot(M456_flare.beta1, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta0 with different sample sizes")

ggplot(M456_flare.beta1, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta0 with different sample sizes")


M456_flare.beta2 <- data.frame(sep=rep(c("M6", "M5", "M4"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M6_100_output$rmse[3], M6_500_output$rmse[3], M6_1000_output$rmse[3],
                                      M5_100_output$rmse[3], M5_500_output$rmse[3], M5_1000_output$rmse[3],
                                      M4_150_output$rmse[3], M4_500_output$rmse[3], M4_1000_output$rmse[3]),
                               bias=c(M6_100_output$bias[3], M6_500_output$bias[3], M6_1000_output$bias[3],
                                      M5_100_output$bias[3], M5_500_output$bias[3], M5_1000_output$bias[3],
                                      M4_150_output$bias[3], M4_500_output$bias[3], M4_1000_output$bias[3]))

ggplot(M456_flare.beta2, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta1 with different sample sizes")

ggplot(M456_flare.beta2, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta1 with different sample sizes")


M456_flare.sigma <- data.frame(sep=rep(c("M6", "M5", "M4"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M6_100_output$rmse[4], M6_500_output$rmse[4], M6_1000_output$rmse[4],
                                      M5_100_output$rmse[4], M5_500_output$rmse[4], M5_1000_output$rmse[4],
                                      M4_150_output$rmse[4], M4_500_output$rmse[4], M4_1000_output$rmse[4]),
                               bias=c(M6_100_output$bias[4], M6_500_output$bias[4], M6_1000_output$bias[4],
                                      M5_100_output$bias[4], M5_500_output$bias[4], M5_1000_output$bias[4],
                                      M4_150_output$bias[4], M4_500_output$bias[4], M4_1000_output$bias[4]))

ggplot(M456_flare.sigma, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of sigma with different sample size")

ggplot(M456_flare.sigma, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of sigma with different sample sizes")


M456_flare.alpha <- data.frame(sep=rep(c("M6", "M5", "M4"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M6_100_output$rmse[5], M6_500_output$rmse[5], M6_1000_output$rmse[5],
                                      M5_100_output$rmse[5], M5_500_output$rmse[5], M5_1000_output$rmse[5],
                                      M4_150_output$rmse[5], M4_500_output$rmse[5], M4_1000_output$rmse[5]),
                               bias=c(M6_100_output$bias[5], M6_500_output$bias[5], M6_1000_output$bias[5],
                                      M5_100_output$bias[5], M5_500_output$bias[5], M5_1000_output$bias[5],
                                      M4_150_output$bias[5], M4_500_output$bias[5], M4_1000_output$bias[5]))

ggplot(M456_flare.alpha, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of alpha with different sample sizes")

ggplot(M456_flare.alpha, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of alpha with different sample sizes")


### Produce Web Table 4 (first 9 rows) and Web Table 5 (first 9 rows)
### Produce Web Figures 7--12 (c), Web Figures 13--14 (a), and Web Figures 15--18 (c)
matout <-  read.table(file = "outM7_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M7_100_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                           sigma = 0.5, alpha = 0.04)


matout <-  read.table(file = "outM7_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M7_500_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                           sigma = 0.5, alpha = 0.04)


matout <-  read.table(file = "outM7_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M7_1000_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.04)


matout <-  read.table(file = "outM8_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M8_100_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                           sigma = 0.5, alpha = 0.2)


matout <-  read.table(file = "outM8_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M8_500_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                           sigma = 0.5, alpha = 0.2)


matout <-  read.table(file = "outM8_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M8_1000_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.2)


matout <-  read.table(file = "outM9_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M9_100_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                           sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM9_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M9_500_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                           sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM9_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M9_1000_output <- rMSE_bias(matout, lambda = c(0.5, 0.5), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.5)


M789_flare.lambda <- data.frame(sep=rep(c("M9", "M8", "M7"), each=3),
                                size=rep(c(100, 500, 1000),3), 
                                rmse=c(M9_100_output$rmse[1], M9_500_output$rmse[1], M9_1000_output$rmse[1],
                                       M8_100_output$rmse[1], M8_500_output$rmse[1], M8_1000_output$rmse[1],
                                       M7_100_output$rmse[1], M7_500_output$rmse[1], M7_1000_output$rmse[1]),
                                bias=c(M9_100_output$bias[1], M9_500_output$bias[1], M9_1000_output$bias[1],
                                       M8_100_output$bias[1], M8_500_output$bias[1], M8_1000_output$bias[1],
                                       M7_100_output$bias[1], M7_500_output$bias[1], M7_1000_output$bias[1]))

ggplot(M789_flare.lambda, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of lambda with different sample size")

ggplot(M789_flare.lambda, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of lambda with different sample sizes")


M789_flare.beta1 <- data.frame(sep=rep(c("M9", "M8", "M7"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M9_100_output$rmse[2], M9_500_output$rmse[2], M9_1000_output$rmse[2],
                                      M8_100_output$rmse[2], M8_500_output$rmse[2], M8_1000_output$rmse[2],
                                      M7_100_output$rmse[2], M7_500_output$rmse[2], M7_1000_output$rmse[2]),
                               bias=c(M9_100_output$bias[2], M9_500_output$bias[2], M9_1000_output$bias[2],
                                      M8_100_output$bias[2], M8_500_output$bias[2], M8_1000_output$bias[2],
                                      M7_100_output$bias[2], M7_500_output$bias[2], M7_1000_output$bias[2]))

ggplot(M789_flare.beta1, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta0 with different sample sizes")

ggplot(M789_flare.beta1, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta0 with different sample sizes")

M789_flare.beta2 <- data.frame(sep=rep(c("M9", "M8", "M7"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M9_100_output$rmse[3], M9_500_output$rmse[3], M9_1000_output$rmse[3],
                                      M8_100_output$rmse[3], M8_500_output$rmse[3], M8_1000_output$rmse[3],
                                      M7_100_output$rmse[3], M7_500_output$rmse[3], M7_1000_output$rmse[3]),
                               bias=c(M9_100_output$bias[3], M9_500_output$bias[3], M9_1000_output$bias[3],
                                      M8_100_output$bias[3], M8_500_output$bias[3], M8_1000_output$bias[3],
                                      M7_100_output$bias[3], M7_500_output$bias[3], M7_1000_output$bias[3]))

ggplot(M789_flare.beta2, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta1 with different sample sizes")

ggplot(M789_flare.beta2, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biasess of beta1 with different sample sizes")

M789_flare.beta3 <- data.frame(sep=rep(c("M9", "M8", "M7"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M9_100_output$rmse[4], M9_500_output$rmse[4], M9_1000_output$rmse[4],
                                      M8_100_output$rmse[4], M8_500_output$rmse[4], M8_1000_output$rmse[4],
                                      M7_100_output$rmse[4], M7_500_output$rmse[4], M7_1000_output$rmse[4]),
                               bias=c(M9_100_output$bias[4], M9_500_output$bias[4], M9_1000_output$bias[4],
                                      M8_100_output$bias[4], M8_500_output$bias[4], M8_1000_output$bias[4],
                                      M7_100_output$bias[4], M7_500_output$bias[4], M7_1000_output$bias[4]))

ggplot(M789_flare.beta3, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta2 with different sample sizes")

ggplot(M789_flare.beta3, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta2 with different sample sizes")


M789_flare.sigma <- data.frame(sep=rep(c("M9", "M8", "M7"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M9_100_output$rmse[5], M9_500_output$rmse[5], M9_1000_output$rmse[5],
                                      M8_100_output$rmse[5], M8_500_output$rmse[5], M8_1000_output$rmse[5],
                                      M7_100_output$rmse[5], M7_500_output$rmse[5], M7_1000_output$rmse[5]),
                               bias=c(M9_100_output$bias[5], M9_500_output$bias[5], M9_1000_output$bias[5],
                                      M8_100_output$bias[5], M8_500_output$bias[5], M8_1000_output$bias[5],
                                      M7_100_output$bias[5], M7_500_output$bias[5], M7_1000_output$bias[5]))

ggplot(M789_flare.sigma, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of sigma with different sample sizes")

ggplot(M789_flare.sigma, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of sigma with different sample sizes")

M789_flare.alpha <- data.frame(sep=rep(c("M9", "M8", "M7"), each=3),
                               size=rep(c(100, 500, 1000),3), 
                               rmse=c(M9_100_output$rmse[6], M9_500_output$rmse[6], M9_1000_output$rmse[6],
                                      M8_100_output$rmse[6], M8_500_output$rmse[6], M8_1000_output$rmse[6],
                                      M7_100_output$rmse[6], M7_500_output$rmse[6], M7_1000_output$rmse[6]),
                               bias=c(M9_100_output$bias[6], M9_500_output$bias[6], M9_1000_output$bias[6],
                                      M8_100_output$bias[6], M8_500_output$bias[6], M8_1000_output$bias[6],
                                      M7_100_output$bias[6], M7_500_output$bias[6], M7_1000_output$bias[6]))

ggplot(M789_flare.alpha, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of alpha with different sample sizes")

ggplot(M789_flare.alpha, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of alpha with different sample sizes")


### Produce Web Table 4 (second 9 rows) and Web Table 5 (second 9 rows)
### Produce Web Figures 7--12 (d), Web Figures 13--14 (b), and Web Figures 15--18 (d)
matout <-  read.table(file = "outM10_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M10_100_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.04)


matout <-  read.table(file = "outM10_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M10_500_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.04)


matout <-  read.table(file = "outM10_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M10_1000_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                             sigma = 0.5, alpha = 0.04)


matout <-  read.table(file = "outM11_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M11_100_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.2)


matout <-  read.table(file = "outM11_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M11_500_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.2)


matout <-  read.table(file = "outM11_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M11_1000_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                             sigma = 0.5, alpha = 0.2)


matout <-  read.table(file = "outM12_100.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M12_100_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM12_500.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M12_500_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                            sigma = 0.5, alpha = 0.5)


matout <-  read.table(file = "outM12_1000.txt")
matout <- as.matrix(matout)
matout <- matout[abs(as.numeric(matout[,"prop1..flare."])-as.numeric(matout[,"prop2..flare."]))<0.99,]

M12_1000_output <- rMSE_bias(matout, lambda = c(0.9, 0.1), beta = c(-2, 1, 13), 
                             sigma = 0.5, alpha = 0.5)


M101112_flare.lambda <- data.frame(sep=rep(c("M12", "M11", "M10"), each=3),
                                   size=rep(c(100, 500, 1000),3), 
                                   rmse=c(M12_100_output$rmse[1], M12_500_output$rmse[1], M12_1000_output$rmse[1],
                                          M11_100_output$rmse[1], M11_500_output$rmse[1], M11_1000_output$rmse[1],
                                          M10_100_output$rmse[1], M10_500_output$rmse[1], M10_1000_output$rmse[1]),
                                   bias=c(M12_100_output$bias[1], M12_500_output$bias[1], M12_1000_output$bias[1],
                                          M11_100_output$bias[1], M11_500_output$bias[1], M11_1000_output$bias[1],
                                          M10_100_output$bias[1], M10_500_output$bias[1], M10_1000_output$bias[1]))

ggplot(M101112_flare.lambda, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of lambda with different sample sizes")

ggplot(M101112_flare.lambda, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of lambda with different sample sizes")


M101112_flare.beta1 <- data.frame(sep=rep(c("M12", "M11", "M10"), each=3),
                                  size=rep(c(100, 500, 1000),3), 
                                  rmse=c(M12_100_output$rmse[2], M12_500_output$rmse[2], M12_1000_output$rmse[2],
                                         M11_100_output$rmse[2], M11_500_output$rmse[2], M11_1000_output$rmse[2],
                                         M10_100_output$rmse[2], M10_500_output$rmse[2], M10_1000_output$rmse[2]),
                                  bias=c(M12_100_output$bias[2], M12_500_output$bias[2], M12_1000_output$bias[2],
                                         M11_100_output$bias[2], M11_500_output$bias[2], M11_1000_output$bias[2],
                                         M10_100_output$bias[2], M10_500_output$bias[2], M10_1000_output$bias[2]))

ggplot(M101112_flare.beta1, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta0 with different sample sizes")

ggplot(M101112_flare.beta1, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta0 with different sample sizes")


M101112_flare.beta2 <- data.frame(sep=rep(c("M12", "M11", "M10"), each=3),
                                  size=rep(c(100, 500, 1000),3), 
                                  rmse=c(M12_100_output$rmse[3], M12_500_output$rmse[3], M12_1000_output$rmse[3],
                                         M11_100_output$rmse[3], M11_500_output$rmse[3], M11_1000_output$rmse[3],
                                         M10_100_output$rmse[3], M10_500_output$rmse[3], M10_1000_output$rmse[3]),
                                  bias=c(M12_100_output$bias[3], M12_500_output$bias[3], M12_1000_output$bias[3],
                                         M11_100_output$bias[3], M11_500_output$bias[3], M11_1000_output$bias[3],
                                         M10_100_output$bias[3], M10_500_output$bias[3], M10_1000_output$bias[3]))

ggplot(M101112_flare.beta2, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta1 with different sample sizes")

ggplot(M101112_flare.beta2, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta1 with different sample sizes")


M101112_flare.beta3 <- data.frame(sep=rep(c("M12", "M11", "M10"), each=3),
                                  size=rep(c(100, 500, 1000),3), 
                                  rmse=c(M12_100_output$rmse[4], M12_500_output$rmse[4], M12_1000_output$rmse[4],
                                         M11_100_output$rmse[4], M11_500_output$rmse[4], M11_1000_output$rmse[4],
                                         M10_100_output$rmse[4], M10_500_output$rmse[4], M10_1000_output$rmse[4]),
                                  bias=c(M12_100_output$bias[4], M12_500_output$bias[4], M12_1000_output$bias[4],
                                         M11_100_output$bias[4], M11_500_output$bias[4], M11_1000_output$bias[4],
                                         M10_100_output$bias[4], M10_500_output$bias[4], M10_1000_output$bias[4]))

ggplot(M101112_flare.beta3, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of beta2 with different sample sizes")

ggplot(M101112_flare.beta3, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of beta2 with different sample sizes")



M101112_flare.sigma <- data.frame(sep=rep(c("M12", "M11", "M10"), each=3),
                                  size=rep(c(100, 500, 1000),3), 
                                  rmse=c(M12_100_output$rmse[5], M12_500_output$rmse[5], M12_1000_output$rmse[5],
                                         M11_100_output$rmse[5], M11_500_output$rmse[5], M11_1000_output$rmse[5],
                                         M10_100_output$rmse[5], M10_500_output$rmse[5], M10_1000_output$rmse[5]),
                                  bias=c(M12_100_output$bias[5], M12_500_output$bias[5], M12_1000_output$bias[5],
                                         M11_100_output$bias[5], M11_500_output$bias[5], M11_1000_output$bias[5],
                                         M10_100_output$bias[5], M10_500_output$bias[5], M10_1000_output$bias[5]))

ggplot(M101112_flare.sigma, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of sigma with different sample size")

ggplot(M101112_flare.sigma, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of sigma with different sample sizes")


M101112_flare.alpha <- data.frame(sep=rep(c("M12", "M11", "M10"), each=3),
                                  size=rep(c(100, 500, 1000),3), 
                                  rmse=c(M12_100_output$rmse[6], M12_500_output$rmse[6], M12_1000_output$rmse[6],
                                         M11_100_output$rmse[6], M11_500_output$rmse[6], M11_1000_output$rmse[6],
                                         M10_100_output$rmse[6], M10_500_output$rmse[6], M10_1000_output$rmse[6]),
                                  bias=c(M12_100_output$bias[6], M12_500_output$bias[6], M12_1000_output$bias[6],
                                         M11_100_output$bias[6], M11_500_output$bias[6], M11_1000_output$bias[6],
                                         M10_100_output$bias[6], M10_500_output$bias[6], M10_1000_output$bias[6]))

ggplot(M101112_flare.alpha, aes(x=size, y=rmse, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("RMSEs of alpha with different sample sizes")

ggplot(M101112_flare.alpha, aes(x=size, y=bias, group=sep)) +
  geom_line(aes(color=sep))+
  geom_point(aes(color=sep))+ 
  ggtitle("Biases of alpha with different sample sizes")


