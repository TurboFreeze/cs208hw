##
## q3regression.r
##
## Differentially private regression
##
## JH 2019/03/10
##



### PART (A)

# optimal upper bound found in previous part
b.optimal <- 15

# data generating process
dgp <- function (n, lambda=10) {
  rpois(n, lambda)
}

# sign function
sgn <- function(x) {
  return(ifelse(x < 0, -1, 1))
}

# clamping helper function
clamp <- function (data, a, b) {
  data.clamped <- data
  data.clamped[data < a] <- a
  data.clamped[data > b] <- b
  data.clamped
}

# laplace random draws
rlap <- function(mu=0, s=1, size=1) {
  p <- runif(size) - 0.5
  draws <- mu - s * sgn(p) * log(1 - 2 * abs(p))
  draws
}

# differentially private mean release
dpmean <- function (data, epsilon, a=0, b) {
  mean.actual <- mean(clamp(data, a, b))
  
  # generate noise by Laplace mechanism
  data.len <- length(data)
  sensitivity <- (b - a) / data.len
  laplace.shift <- sensitivity / epsilon
  noise <- rlap(s=laplace.shift, size=1)
  
  # inject noise
  mean.noisy <- mean.actual + noise
  
  # apply clamping
  mean.clamped <- clamp(mean.noisy, a, b)
}


# Differentially private regression slope release (from provided code)
dpslope <- function (y, x, ylower=0, yupper=b.optimal, xlower=0, xupper=b.optimal, eps1, eps2){
  x <- clamp(x, xlower, xupper)
  y <- clamp(y, ylower, yupper)
  
  n <- length(x)
  sens.sxy <- (xupper - xlower) * (yupper - ylower)
  sens.sxx  <- (xupper - xlower)^2
  
  scale.sxy <- sens.sxy / eps1
  scale.sxx <- sens.sxx / eps2
  
  sensitive.value <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2) 
  
  release.sxy <- sum((x - mean(x)) * (y - mean(y))) + rlap(mu=0, s=scale.sxy, size=1)
  release.sxx <- sum((x - mean(x))^2) + rlap(mu=0, s=scale.sxx, size=1)
  
  postprocess.beta <- release.sxy / release.sxx
  postprocess.beta
}

# Differentially private linear regression release of intercept and slope
dpregression <- function (y, x, ylower=0, yupper=b.optimal, xlower=0, xupper=b.optimal,
                          epsilon, epsilon.partition=c(0.25, 0.25, 0.25, 0.25)) {
  x <- clamp(x, xlower, xupper)
  y <- clamp(y, ylower, yupper)
  
  # calculate dp slope
  dpbeta <- dpslope(y, x, ylower, yupper, xlower, xupper,
                    epsilon * epsilon.partition[1], epsilon * epsilon.partition[2])
  
  # calculate dp intercept
  dpmean.x <- dpmean(x, epsilon * epsilon.partition[3], b=xupper)
  dpmean.y <- dpmean(y, epsilon * epsilon.partition[4], b=yupper)
  dpalpha <- dpmean.y - dpbeta * dpmean.x
  
  list(dpalpha=dpalpha, dpbeta=dpbeta)
}



### PART (B)

n <- 1000
alpha <- 1
beta <- 1
sigma <- 1
epsilon <- 1

# mean squared residuals for quantifying utility
msr <- function (y, x, beta.hat, alpha.hat) {
  sum((y - beta.hat * x - alpha.hat)^2) / length(y)
}

simulate <- function(n.trials=100) {
  msrs.dp <- vector("numeric", n.trials)
  msrs.non <- vector("numeric", n.trials)
  for (i in 1:n.trials) {
    # generate data
    x <- dgp(n, lambda=10)
    y <- beta * x + alpha + rnorm(n, 0, sigma)
    
    # calculate differentially private release
    results.dp <- dpregression(y, x, epsilon=epsilon)
    msrs.dp[i] <- msr(y, x, results.dp$dpbeta, results.dp$dpalpha)
    
    # calculate non-private release
    results.non <- lm(y ~ x)
    msrs.non[i] <- msr(y, x, coef(results.non)[[2]], coef(results.non)[[1]])
  }
  data.frame(n=(1:n.trials), dp=msrs.dp, non=msrs.non)
}

msrs.data <- simulate(1000)

# visualize
library(ggplot2)
y.window <- c(0, 10)

# monte carlo scatterplot
msrs.scatterplot <- ggplot(msrs.data, aes(x=n)) +
  geom_point(aes(y=dp), color="red", alpha=0.7, size=0.5) +
  geom_point(aes(y=non), color="green", alpha=0.7, size=0.5) +
  ylim(y.window) +
  ylab("mean squared residuals") +
  theme_bw(); msrs.scatterplot

# monte carlo traceplot
msrs.traceplot <- ggplot(msrs.data, aes(x=n)) +
  geom_line(aes(y=dp), color="red") +
  geom_line(aes(y=non), color="green") +
  ylim(y.window) +
  ylab("mean squared residuals") +
  theme_bw(); msrs.traceplot

# distributions
msrs.distribution <- ggplot(msrs.data) +
  xlim(c(0, 5)) +
  xlab("mean squared residuals") +
  geom_density(aes(dp), fill="red", alpha=0.7) +
  geom_density(aes(non), fill="green", alpha=0.7) +
  theme_bw(); msrs.distribution
# distribution histograms
#hist(msrs.data$dp[msrs.data$dp < 10]) # uncomment if no ggplot

library(gridExtra)
gridplots <- grid.arrange(msrs.scatterplot, msrs.traceplot)

ggsave("msrs.jpg", gridplots)
ggsave("msrsdist.jpg", msrs.distribution)



### PART (C)
# helper function to normalize to distribution of weights
normalize <- function(vec) {
  vec / sum(vec)
}

# hyperparameter space
param.domain <- seq(0.1, 1, by=0.1)
# recursive function to generate all hyperparamater combinations
param.combs <- function (hypers) {
  if (length(hypers) == 4) {
    normalize(hypers)
  } else {
    res <- c()
    for (i in 1:length(param.domain)) {
      res <- rbind(res, param.combs(c(hypers, param.domain[i])))
    }
    res
  }
}
# generate hyperparameter combinations
param.space <- param.combs(c())
param.msr <- vector("numeric", nrow(param.space))

# calculate mean squared residuals
for (i in 1:nrow(param.space)) {
  # generate data
  x <- dgp(n, lambda=10)
  y <- beta * x + alpha + rnorm(n, 0, sigma)
  
  # calculate differentially private release
  results.dp <- dpregression(y, x, epsilon=epsilon, epsilon.partition=param.space[i,])
  param.msr[i] <- msr(y, x, results.dp$dpbeta, results.dp$dpalpha)
}

# equal partition of epsilon
results.dp <- dpregression(y, x, epsilon=epsilon, epsilon.partition=param.space[i,])
equal.msr <- msr(y, x, results.dp$dpbeta, results.dp$dpalpha); equal.msr

# get minimum mean squared residuals
min.index <- c()
min.space <- c()
min.msr <- c()
param.msr.working <- param.msr
param.space.working <- param.space
# get top 10 hyperparameter configurations
for (i in 1:10) {
  m <- which.min(param.msr.working)
  min.index <- c(min.index, m)
  min.msr <- c(min.msr, param.msr.working[m])
  min.space <- rbind(min.space, param.space.working[m,])
  param.msr.working <- param.msr.working[-1 * m]
  param.space.working <- param.space.working[-1 * m,]
}
# find mean and median of hyperparameters
min.msr
mea <- apply(min.space, mean, MARGIN = 2); mea
med <- apply(min.space, median, MARGIN = 2); med
