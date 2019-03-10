##
## q2.r
##
## Evaluating DP algorithms with synthetic data
##
## JH 2019/03/10
##



##3 PART (A)
# data generating process
dgp <- function (n, lambda=10) {
  rpois(n, lambda)
}



### PART (B)

# sign function
sgn <- function(x) {
  return(ifelse(x < 0, -1, 1))
}

# laplace random draws
rlap <- function(mu=0, s=1, size=1) {
  p <- runif(size) - 0.5
  draws <- mu - s * sgn(p) * log(1 - 2 * abs(p))
  draws
}

# clamping helper function
clamp <- function (data, a, b) {
  data.clamped <- data
  data.clamped[data < a] <- a
  data.clamped[data > b] <- b
  data.clamped
}

# differentially private mechanism (clamping)
dpclamping <- function (data, epsilon, a=0, b) {
  mean.actual <- mean(data)
  
  # generate noise by Laplace mechanism
  data.len <- length(data)
  laplace.shift <- 2 / data.len
  noise <- rlap(s=laplace.shift, size=1)
  
  # inject noise
  mean.noisy <- mean.actual + noise
  
  # apply clamping
  mean.clamped <- clamp(mean.noisy, a, b)
}



### PART (C)

# parameters
n <- 200
epsilon <- 0.5
b.seq <- seq(0, 30, by=0.1)
b.num <- length(b.seq)

# generate data
data.poisson <- dgp(n, lambda=10)

# calculate RMSE
rmse <- function (data, m) {
  sqrt(sum((data - m)^2))
}

rmses <- vector("numeric", b.num)
for (i in 1:b.num) {
  # calculate dp query of mean
  dpmean <- dpclamping(data=data.poisson, epsilon=epsilon, b=b.seq[i])
  # calculate and store RMSE
  rmses[i] <- rmse(data.poisson, dpmean)
}

# find index of minimum RMSE (optimal value of b)
b.optimal <- b.seq[which.min(rmses)]

#plot(b.seq, rmses, type='l') # uncomment if no ggplot

# visualize
rmsedata <- data.frame(b=b.seq, rmse=rmses)
library(ggplot2)
clamping.plot <- ggplot(rmsedata, aes(x=b, y=rmse)) +
  geom_line() +
  geom_vline(xintercept=b.optimal, color="red", alpha=0.8, linetype="dashed") +
  theme_bw()
ggsave("clamping.jpg", clamping.plot)
