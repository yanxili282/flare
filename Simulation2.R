#### Author: Yanxi Li, yli282@g.uky.edu
### Required packages
library(MASS)
library(emg)
library(mixtools)
library(ggplot2)
library(plotly)
library("gg3D")
library(scatterplot3d)
library(plotly)

### load functions:
sim.flare_reg2 <- function(xmin, xmax, n, beta, lambda, sigma, alpha){
  ## x is a matrix without intercept
  m = length(beta)
  X <- matrix(runif(n*(m-1), min = xmin, max = xmax), 
              nrow = n, ncol = (m-1), byrow = FALSE, dimnames = NULL)
  X <- cbind(1,X)
  lambda1 = lambda[1]
  
  E <- c()
  G <- c()
  for (i in 1:n){
    indicator <- runif(1, min = 0, max = 1)
    if (indicator <= lambda1){
      eps <- rnorm(1, mean=0, sd=sigma)
      g <- "1"
      E <- c(E, eps)
      G <- c(G, g)
    } else{
      eps <- rexp(1, rate=alpha)
      g <- "2"
      E <- c(E, eps)
      G <- c(G, g)
    }
  }
  
  E <- matrix(data = E, nrow = n, ncol = 1, byrow = FALSE,
              dimnames = NULL)
  
  Y <- X %*% beta + E
  
  X <- X[,-1]
  
  sim.out <- list("X" = X,
                  "Y" = Y,
                  "Group"=G)
  
  return(sim.out)
}

### Produce Web Figure 3 in the Supporting Information:
## M1_100:
out_M5_100 <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                             lambda=c(1/3, 2/3), sigma=1/2, alpha=1/20)

df_M5_100 <- data.frame(x = out_M5_100$X, y = out_M5_100$Y, group = out_M5_100$Group)

ggplot(df_M5_100, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Well-separated", subtitle = "n=100")

## M1_500:
out_M5_500 <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                             lambda=c(1/3, 2/3), sigma=1/2, alpha=1/20)

df_M5_500 <- data.frame(x = out_M5_500$X, y = out_M5_500$Y, group = out_M5_500$Group)

ggplot(df_M5_500, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Well-separated", subtitle = "n=500")

## M1_1000:
out_M5_1000 <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                              lambda=c(1/3, 2/3), sigma=1/2, alpha=1/20)

df_M5_1000 <- data.frame(x = out_M5_1000$X, y = out_M5_1000$Y, group = out_M5_1000$Group)

ggplot(df_M5_1000, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Well-separated", subtitle = "n=1000")

## M2_100:
out_M5_100_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                    lambda=c(1/3, 2/3), sigma=1/2, alpha=1/6)

df_M5_100_medium <- data.frame(x = out_M5_100_medium$X, y = out_M5_100_medium$Y, group = out_M5_100_medium$Group)

ggplot(df_M5_100_medium, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Moderately-separated", subtitle = "n=100")

## M2_500:
out_M5_500_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                    lambda=c(1/3, 2/3), sigma=1/2, alpha=1/6)

df_M5_500_medium <- data.frame(x = out_M5_500_medium$X, y = out_M5_500_medium$Y, group = out_M5_500_medium$Group)

ggplot(df_M5_500_medium, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Moderately-separated", subtitle = "n=500")

## M2_1000:
out_M5_1000_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                     lambda=c(1/3, 2/3), sigma=1/2, alpha=1/6)

df_M5_1000_medium <- data.frame(x = out_M5_1000_medium$X, y = out_M5_1000_medium$Y, group = out_M5_1000_medium$Group)

ggplot(df_M5_1000_medium, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Moderately-separated", subtitle = "n=1000")

## M3_100:
out_M5_100_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                  lambda=c(1/3, 2/3), sigma=1/2, alpha=1/2)

df_M5_100_poor <- data.frame(x = out_M5_100_poor$X, y = out_M5_100_poor$Y, group = out_M5_100_poor$Group)

ggplot(df_M5_100_poor, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Overlapping", subtitle = "n=100")

## M3_500:
out_M5_500_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                  lambda=c(1/3, 2/3), sigma=1/2, alpha=1/2)

df_M5_500_poor <- data.frame(x = out_M5_500_poor$X, y = out_M5_500_poor$Y, group = out_M5_500_poor$Group)

ggplot(df_M5_500_poor, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Overlapping", subtitle = "n=500")

## M3_1000:
out_M5_1000_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                   lambda=c(1/3, 2/3), sigma=1/2, alpha=1/2)

df_M5_1000_poor <- data.frame(x = out_M5_1000_poor$X, y = out_M5_1000_poor$Y, group = out_M5_1000_poor$Group)

ggplot(df_M5_1000_poor, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Overlapping", subtitle = "n=1000")

### Produce Web Figure 4 in the Supporting Information:
## M4_100:
out_M7_100 <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                             lambda=c(0.9, 0.1), sigma=1/2, alpha=1/20)

df_M7_100 <- data.frame(x = out_M7_100$X, y = out_M7_100$Y, group = out_M7_100$Group)

ggplot(df_M7_100, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Well-separated", subtitle = "n=100")

## M4_500:
out_M7_500 <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                             lambda=c(0.9, 0.1), sigma=1/2, alpha=1/20)

df_M7_500 <- data.frame(x = out_M7_500$X, y = out_M7_500$Y, group = out_M7_500$Group)

ggplot(df_M7_500, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Well-separated", subtitle = "n=500")

## M4_1000:
out_M7_1000 <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                              lambda=c(0.9, 0.1), sigma=1/2, alpha=1/20)

df_M7_1000 <- data.frame(x = out_M7_1000$X, y = out_M7_1000$Y, group = out_M7_1000$Group)

ggplot(df_M7_1000, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Well-separated", subtitle = "n=1000")


## M5_100:
out_M7_100_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                    lambda=c(0.9, 0.1), sigma=1/2, alpha=1/6)

df_M7_100_medium <- data.frame(x = out_M7_100_medium$X, y = out_M7_100_medium$Y, group = out_M7_100_medium$Group)

ggplot(df_M7_100_medium, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Moderately-separated", subtitle = "n=100")

## M5_500:
out_M7_500_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                    lambda=c(0.9, 0.1), sigma=1/2, alpha=1/6)

df_M7_500_medium <- data.frame(x = out_M7_500_medium$X, y = out_M7_500_medium$Y, group = out_M7_500_medium$Group)

ggplot(df_M7_500_medium, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Moderately-separated", subtitle = "n=500")

## M5_1000:
out_M7_1000_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                     lambda=c(0.9, 0.1), sigma=1/2, alpha=1/6)

df_M7_1000_medium <- data.frame(x = out_M7_1000_medium$X, y = out_M7_1000_medium$Y, group = out_M7_1000_medium$Group)

ggplot(df_M7_1000_medium, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Moderately-separated", subtitle = "n=1000")

## M6_100:
out_M7_100_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(9, 3), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)

df_M7_100_poor <- data.frame(x = out_M7_100_poor$X, y = out_M7_100_poor$Y, group = out_M7_100_poor$Group)

ggplot(df_M7_100_poor, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Overlapping", subtitle = "n=100")

## M6_500:
out_M7_500_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(9, 3), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)

df_M7_500_poor <- data.frame(x = out_M7_500_poor$X, y = out_M7_500_poor$Y, group = out_M7_500_poor$Group)

ggplot(df_M7_500_poor, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Overlapping", subtitle = "n=500")

## M6_1000:
out_M7_1000_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(9, 3), 
                                   lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)

df_M7_1000_poor <- data.frame(x = out_M7_1000_poor$X, y = out_M7_1000_poor$Y, group = out_M7_1000_poor$Group)

ggplot(df_M7_1000_poor, aes(x = x, y = y, color = group)) +
  geom_point()+ 
  labs(title = "Overlapping", subtitle = "n=1000")

### Produce Web Figure 5 in the Supporting Information:
## M7_100:
out_M6_100 <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                             lambda=c(1/2, 1/2), sigma=1/2, alpha=1/25)

df_M6_100 <- data.frame(x = out_M6_100$X[,1], y = out_M6_100$X[,2], z=out_M6_100$Y, group = out_M6_100$Group)

colors <- c("#E69F00", "#960018")
colors <- colors[as.numeric(df_M6_100$group)]

plot_ly(df_M6_100 , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Well-separated (n=100)",
                                                                      scene = list(xaxis = list(title = 'x1'),
                                                                                   yaxis = list(title = 'x2'),
                                                                                   zaxis = list(title = 'y')))

## M7_500:
out_M6_500 <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                             lambda=c(1/2, 1/2), sigma=1/2, alpha=1/25)

df_M6_500 <- data.frame(x = out_M6_500$X[,1], y = out_M6_500$X[,2], z=out_M6_500$Y, group = out_M6_500$Group)

colors <- c( "#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_500$group)]
plot_ly(df_M6_500 , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Well-separated (n=500)",
                                                                      scene = list(xaxis = list(title = 'x1'),
                                                                                   yaxis = list(title = 'x2'),
                                                                                   zaxis = list(title = 'y')))

## M7_1000:
out_M6_1000 <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                              lambda=c(1/2, 1/2), sigma=1/2, alpha=1/25)

df_M6_1000 <- data.frame(x = out_M6_1000$X[,1], y = out_M6_1000$X[,2], z=out_M6_1000$Y, group = out_M6_1000$Group)

colors <- c( "#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_1000$group)]

plot_ly(df_M6_1000 , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Well-separated (n=1000)",
                                                                       scene = list(xaxis = list(title = 'x1'),
                                                                                    yaxis = list(title = 'x2'),
                                                                                    zaxis = list(title = 'y')))

## M8_100:
out_M6_100_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                    lambda=c(1/2, 1/2), sigma=1/2, alpha=0.2)

df_M6_100_medium <- data.frame(x = out_M6_100_medium$X[,1], y = out_M6_100_medium$X[,2], z=out_M6_100_medium$Y, group = out_M6_100_medium$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_100_medium$group)]

plot_ly(df_M6_100_medium , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Moderately-separated (n=100)",
                                                                             scene = list(xaxis = list(title = 'x1'),
                                                                                          yaxis = list(title = 'x2'),
                                                                                          zaxis = list(title = 'y')))

## M8_500:
out_M6_500_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                    lambda=c(1/2, 1/2), sigma=1/2, alpha=0.2)

df_M6_500_medium <- data.frame(x = out_M6_500_medium$X[,1], y = out_M6_500_medium$X[,2], z=out_M6_500_medium$Y, group = out_M6_500_medium$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_500_medium$group)]

plot_ly(df_M6_500_medium , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Moderately-separated (n=500)",
                                                                             scene = list(xaxis = list(title = 'x1'),
                                                                                          yaxis = list(title = 'x2'),
                                                                                          zaxis = list(title = 'y')))


## M8_1000:
out_M6_1000_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                     lambda=c(1/2, 1/2), sigma=1/2, alpha=0.2)

df_M6_1000_medium <- data.frame(x = out_M6_1000_medium$X[,1], y = out_M6_1000_medium$X[,2], z=out_M6_1000_medium$Y, group = out_M6_1000_medium$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_1000_medium$group)]

plot_ly(df_M6_1000_medium , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Moderately-separated (n=1000)",
                                                                              scene = list(xaxis = list(title = 'x1'),
                                                                                           yaxis = list(title = 'x2'),
                                                                                           zaxis = list(title = 'y')))
## M9_100:
out_M6_100_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                  lambda=c(1/2, 1/2), sigma=1/2, alpha=1/2)

df_M6_100_poor <- data.frame(x = out_M6_100_poor$X[,1], y = out_M6_100_poor$X[,2], z=out_M6_100_poor$Y, group = out_M6_100_poor$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_100_poor$group)]

plot_ly(df_M6_100_poor , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Overlapping (n=100)",
                                                                           scene = list(xaxis = list(title = 'x1'),
                                                                                        yaxis = list(title = 'x2'),
                                                                                        zaxis = list(title = 'y')))


## M9_500:
out_M6_500_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                  lambda=c(1/2, 1/2), sigma=1/2, alpha=1/2)

df_M6_500_poor <- data.frame(x = out_M6_500_poor$X[,1], y = out_M6_500_poor$X[,2], z=out_M6_500_poor$Y, group = out_M6_500_poor$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_500_poor$group)]

plot_ly(df_M6_500_poor , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Overlapping (n=500)",
                                                                           scene = list(xaxis = list(title = 'x1'),
                                                                                        yaxis = list(title = 'x2'),
                                                                                        zaxis = list(title = 'y')))

## M9_1000:
out_M6_1000_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                   lambda=c(1/2, 1/2), sigma=1/2, alpha=1/2)

df_M6_1000_poor <- data.frame(x = out_M6_1000_poor$X[,1], y = out_M6_1000_poor$X[,2], z=out_M6_1000_poor$Y, group = out_M6_1000_poor$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M6_1000_poor$group)]

plot_ly(df_M6_1000_poor , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Overlapping (n=1000)",
                                                                            scene = list(xaxis = list(title = 'x1'),
                                                                                         yaxis = list(title = 'x2'),
                                                                                         zaxis = list(title = 'y')))

### Produce Web Figure 6 in the Supporting Information:
## M10_100:
out_M8_100 <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                             lambda=c(0.9, 0.1), sigma=1/2, alpha=1/25)

df_M8_100 <- data.frame(x = out_M8_100$X[,1], y = out_M8_100$X[,2], z=out_M8_100$Y, group = out_M8_100$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_100$group)]

plot_ly(df_M8_100 , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Well-separated (n=100)",
                                                                      scene = list(xaxis = list(title = 'x1'),
                                                                                   yaxis = list(title = 'x2'),
                                                                                   zaxis = list(title = 'y')))


## M10_500:
out_M8_500 <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                             lambda=c(0.9, 0.1), sigma=1/2, alpha=1/25)

df_M8_500 <- data.frame(x = out_M8_500$X[,1], y = out_M8_500$X[,2], z=out_M8_500$Y, group = out_M8_500$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_500$group)]

plot_ly(df_M8_500 , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Well-separated (n=500)",
                                                                      scene = list(xaxis = list(title = 'x1'),
                                                                                   yaxis = list(title = 'x2'),
                                                                                   zaxis = list(title = 'y')))



## M10_1000:
out_M8_1000 <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                              lambda=c(0.9, 0.1), sigma=1/2, alpha=1/25)

df_M8_1000 <- data.frame(x = out_M8_1000$X[,1], y = out_M8_1000$X[,2], z=out_M8_1000$Y, group = out_M8_1000$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_1000$group)]

plot_ly(df_M8_1000 , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Well-separated (n=1000)",
                                                                       scene = list(xaxis = list(title = 'x1'),
                                                                                    yaxis = list(title = 'x2'),
                                                                                    zaxis = list(title = 'y')))

## M11_100:
out_M8_100_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                    lambda=c(0.9, 0.1), sigma=1/2, alpha=0.2)

df_M8_100_medium <- data.frame(x = out_M8_100_medium$X[,1], y = out_M8_100_medium$X[,2], z=out_M8_100_medium$Y, group = out_M8_100_medium$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_100_medium$group)]

plot_ly(df_M8_100_medium , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Moderately-separated (n=100)",
                                                                             scene = list(xaxis = list(title = 'x1'),
                                                                                          yaxis = list(title = 'x2'),
                                                                                          zaxis = list(title = 'y')))

## M11_500:
out_M8_500_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                    lambda=c(0.9, 0.1), sigma=1/2, alpha=0.2)

df_M8_500_medium <- data.frame(x = out_M8_500_medium$X[,1], y = out_M8_500_medium$X[,2], z=out_M8_500_medium$Y, group = out_M8_500_medium$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_500_medium$group)]

plot_ly(df_M8_500_medium , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Moderately-separated (n=500)",
                                                                             scene = list(xaxis = list(title = 'x1'),
                                                                                          yaxis = list(title = 'x2'),
                                                                                          zaxis = list(title = 'y')))


## M11_1000:
out_M8_1000_medium <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                     lambda=c(0.9, 0.1), sigma=1/2, alpha=0.2)

df_M8_1000_medium <- data.frame(x = out_M8_1000_medium$X[,1], y = out_M8_1000_medium$X[,2], z=out_M8_1000_medium$Y, group = out_M8_1000_medium$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_1000_medium$group)]

plot_ly(df_M8_1000_medium , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Moderately-separated (n=1000)",
                                                                              scene = list(xaxis = list(title = 'x1'),
                                                                                           yaxis = list(title = 'x2'),
                                                                                           zaxis = list(title = 'y')))


## M12_100:
out_M8_100_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=100, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)

df_M8_100_poor <- data.frame(x = out_M8_100_poor$X[,1], y = out_M8_100_poor$X[,2], z=out_M8_100_poor$Y, group = out_M8_100_poor$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_100_poor$group)]

plot_ly(df_M8_100_poor , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Overlapping (n=100)",
                                                                           scene = list(xaxis = list(title = 'x1'),
                                                                                        yaxis = list(title = 'x2'),
                                                                                        zaxis = list(title = 'y')))



## M12_500:
out_M8_500_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=500, beta=c(-2, 1, 13), 
                                  lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)

df_M8_500_poor <- data.frame(x = out_M8_500_poor$X[,1], y = out_M8_500_poor$X[,2], z=out_M8_500_poor$Y, group = out_M8_500_poor$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_500_poor$group)]

plot_ly(df_M8_500_poor , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Overlapping (n=500)",
                                                                           scene = list(xaxis = list(title = 'x1'),
                                                                                        yaxis = list(title = 'x2'),
                                                                                        zaxis = list(title = 'y')))



## M12_1000:
out_M8_1000_poor <- sim.flare_reg2(xmin=-10, xmax=10, n=1000, beta=c(-2, 1, 13), 
                                   lambda=c(0.9, 0.1), sigma=1/2, alpha=1/2)

df_M8_1000_poor <- data.frame(x = out_M8_1000_poor$X[,1], y = out_M8_1000_poor$X[,2], z=out_M8_1000_poor$Y, group = out_M8_1000_poor$Group)

colors <- c("#E69F00", "#56B4E9")
colors <- colors[as.numeric(df_M8_1000_poor$group)]

plot_ly(df_M8_1000_poor , x = ~x, y = ~y, z = ~z, color = ~group)%>% layout(title="Overlapping (n=1000)",
                                                                            scene = list(xaxis = list(title = 'x1'),
                                                                                         yaxis = list(title = 'x2'),
                                                                                         zaxis = list(title = 'y')))


