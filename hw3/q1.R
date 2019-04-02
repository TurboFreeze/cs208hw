## cs208 hw3
## q1.r
##
## Tails, trimming, and Winsorization
##
## JH 2019/03/19
##


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

# clipping helper function
clip <- function(x, lower, upper){
  x.clipped <- x
  x.clipped[x.clipped < lower] <- lower
  x.clipped[x.clipped > upper] <- upper
  x.clipped
}

# helper function for calculating RMSE
rmse <- function (actual, release) {
  sqrt(sum((actual - release) ^ 2))
}

### PART (A)
# trimmed mean
release.mean.trimmed <- function (data, D, epsilon) {
  n <- length(data)
  data.ranks <- rank(data)
  data.trimmed <- data[data.ranks <= (0.95 * n) & data.ranks >= (0.05 * n)]
  m <- (1 / (0.9 * n)) * sum(data.trimmed) + rlap(0, D / (epsilon * n), 1)
  m
}


### PART (C)

# exponential percentile
release.percentile <- function (data, t, lower, upper, nbins=0, epsilon) {
  n <- length(data)
  percentile <- floor(t * n / 100)
  # binning data
  bins <- seq(lower, upper, ceiling((upper - lower) / nbins) + 1)
  if (nbins == 0) {
    bins <- floor(lower):ceiling(upper)
    nbins <- length(bins)
  }
  data.clipped <- clip(data, lower, upper)
  
  # calculate utility function values for each bin
  utility <- rep(NA, nbins)
  for(i in 1:length(utility)){
    utility[i] <- min(-abs(sum(data.clipped <= bins[i]) - percentile),
                          -abs(sum(data.clipped >= bins[i]) - percentile))
  }
  likelihoods <- exp(epsilon * utility / 2)
  probabilities <- likelihoods / sum(likelihoods)
  
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities)
  dprelease <- min(bins[flag]) 
  
  dprelease
}


### PART (D)
# trimmed mean with percentile estimates
release.mean.percestimates <- function (data, D, epsilon, nbins=1000) {
  n <- length(data)
  # dp percentile estimates
  percentile.05 <- release.percentile(data, 5, 0, D, nbins, epsilon / 3)
  percentile.95 <- release.percentile(data, 95, 0, D, nbins, epsilon / 3)
  # calculate dp mean based on these estimates
  data.trimmed <- clip(data, percentile.05, percentile.95)
  lscale <- 3 * (percentile.95 - percentile.05) / (0.9 * epsilon * n)
  m <- (1 / (0.9 * n)) * sum(data.trimmed) + rlap(0, lscale, 1)
  m
}



### PART (F)
# PUMS dataset
data.pums <- read.csv('MaPUMS5full.csv')
epsilon <- 1
D <- 1000000

# ordinary DP mean release
release.mean <- function(data, D, epsilon) {
  lower <- 0
  upper <- D
  n <- length(data)
  
  # calculating scale for Laplace noise
  sensitivity <- (upper - lower) / n
  scale <- sensitivity / epsilon
  
  # DP mean release
  data.clipped <- clip(data, lower, upper)
  sensitive.value <- mean(data.clipped)
  m <- sensitive.value + rlap(mu=0, s=scale, size=1)
  
  m
}

clamped.mean <- function(data, D) {
  clip(mean(data), 0, D)
}

# explore pumas
pumas <- unique(data.pums$puma)

# calculate DP means for each puma
simulate.puma <- function (data, algorithm, n.sims) {
  means.dp <- vector("numeric", n.sims)
  # run simulations and store RMSE
  for (i in 1:n.sims) {
    means.dp[i] <- algorithm(data, D, epsilon)
  }
  means.dp
}

n.sims <- 100
results <- data.frame(puma=character(), algo=character(), dp=numeric(),
                      actual=numeric(), rmse=numeric())
# simulate for each puma
for (puma.id in pumas) {
  # extract data and simulate
  puma.data <- data.pums$income[data.pums$puma == puma.id]
  mean.actual <- clamped.mean(puma.data, D)
  
  # ordinary Laplace mechanism
  mean.ordinary <- simulate.puma(puma.data, release.mean, n.sims)
  rmse.ordinary <- rmse(mean.ordinary, mean.actual)
  results.ordinary <- data.frame(puma=rep(puma.id, n.sims), algo=rep("ordinary", n.sims),
                                dp=mean.ordinary, actual=rep(mean.actual, n.sims),
                                rmse=rmse.ordinary)
  
  # with percentile estimates
  mean.perc <- simulate.puma(puma.data, release.mean.percestimates, n.sims)
  rmse.perc <- rmse(mean.perc, mean.actual)
  results.perc <- data.frame(puma=rep(puma.id, n.sims), algo=rep("perc", n.sims),
                                 dp=mean.perc, actual=rep(mean.actual, n.sims),
                                 rmse=rmse.perc)
  results <- rbind(results, results.ordinary, results.perc)
}

# plotting
library(ggplot2)

# plotting rmse
plot.rmse <- ggplot(results, aes(x=reorder(factor(puma), actual))) +
  geom_point(aes(y=rmse, color=algo)) +
  ggtitle("RMSEs for two DP-mean mechanisms by PUMA") +
  scale_x_discrete(breaks = levels(reorder(factor(results$puma), results$actual))[c(TRUE, rep(FALSE, 2))]) +
  xlab("puma") +
  theme_bw()
plot.rmse
ggsave("rmse.jpg", plot.rmse)


# box and whisker plots
plot.box <- ggplot(results, aes(x=reorder(factor(puma), actual))) + 
  geom_boxplot(aes(y=dp, color=algo), outlier.shape = NA) +
  geom_point(aes(y=actual, group="actual"), color="blue", alpha=0.5, size=0.5) +
  ggtitle("Box and Whisker of DP-means") +
  scale_x_discrete(breaks = levels(reorder(factor(results$puma), results$actual))[c(TRUE, rep(FALSE, 2))]) +
  xlab("puma") +
  ylab("dp mean") +
  theme_bw()
plot.box
ggsave("box.jpg", plot.box)

