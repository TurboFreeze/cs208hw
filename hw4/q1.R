## cs208 hw4
## q1.r
##
## Learning conjunctions in SQ model
##
## JH 2019/04/15
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

# centralized DP mean release
release.mean <- function(data, lower, upper, epsilon) {
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

## Locally differentially private randomized response
release.local <- function(x, values=c(-1,1), epsilon){
  draw <- runif(n=1, min=0, max=1)
  cutoff <- 1 / (1 + exp(epsilon))
  if(draw < cutoff){
    values[values != x]
  }else{
    x
  }
}



### PART (A) implementing centralized and local DP SQ algorithm
sq.centralized <- function (X, y, t, epsilon) {
  d <- ncol(X)
  n <- nrow(X)
  # check predicate x[j] = 0 and y = 1
  X.modified <- abs(X - 1) & y
  # calculate centralized dp means
  prob.estimates <- apply(X.modified, MARGIN = 2, release.mean, 0, 1, epsilon / d) / n
  # output subset estimate based on threshold
  result <- 1:d[prob.estimates <= t]
  result
}

sq.local <- function (X, y, t, epsilon) {
  d <- ncol(X)
  n <- nrow(X)
  # check predicate x[j] = 0 and y = 1
  X.modified <- abs(X - 1) & y
  X.modified.estimates <- apply(X.modified, release.local, values=c(0, 1), epsilon / (n * d))
  # calculate local dp means
  prob.estimates <- apply(X.modified, MARGIN = 2, mean)
  # output subset estimate based on threshold
  result <- 1:d[prob.estimates <= t]
  result
}

### PART (C) classifier learning with DP SQ algorithm

data.pums <- read.csv('CaPUMS5Full.csv')
