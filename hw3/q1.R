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

### PART (A)
# trimmed mean
dpmean.trimmed <- function (n, data, D, epsilon) {
  data.ranks <- rank(data)
  data.trimmed <- data[data.ranks <= (0.95 * n) & data.ranks >= (0.05 * n)]
  m <- (1 / (0.9 * n)) * sum(data.trimmed) + rlap(D / (epsilon * n))
  m
}


### PART (C)

# exponential percentile
release.percentile <- function (data, t, lower, upper, nbins=0, epsilon) {
  n <- length(data)
  percentile <- floor(t * n / 100)
  # binning data
  if (nbins == 0) {
    bins <- floor(lower):ceiling(upper)
    nbins <- length(bins)
  }
  data.clipped <- clip(data, lower, upper)
  # actual percentile value
  sensitive.value <- data[rank(data) == percentile]
  
  # calculate utility function values for each bin
  utility <- rep(NA, nbins)
  for(i in 1:length(utility)){
    utility[i] <- n - max(-abs(sum(data.clipped <= bins[i]) - percentile),
                          -abs(sum(x.clipped>=bins[i]) - percentile))
  }
  likelihoods <- exp(epsilon * utility / 2)
  probabilities <- likelihoods / sum(likelihoods)
  
  flag <- runif(n=1, min=0, max=1) < cumsum(probabilities)
  dprelease <- min(bins[flag]) 
  
  dprelease
}


### PART (D)
# trimmed mean



### PART (F)
# PUMS dataset
data.pums <- read.csv('MaPUMS5full.csv')
epsilon <- 1
D <- 1000000