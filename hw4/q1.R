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
  prob.estimates <- apply(X.modified, MARGIN = 2, release.mean, 0, 1, epsilon / d)
  # output subset estimate based on threshold
  result <- (1:d)[prob.estimates <= t]
  result
}

sq.local <- function (X, y, t, epsilon) {
  d <- ncol(X)
  n <- nrow(X)
  # check predicate x[j] = 0 and y = 1
  X.modified <- abs(X - 1) & y
  X.modified.estimates <- apply(X.modified, MARGIN = c(1, 2), release.local, values=c(0, 1), epsilon=epsilon / (d))
  # calculate local dp means
  prob.estimates <- apply(X.modified.estimates, MARGIN = 2, mean)
  # output subset estimate based on threshold
  result <- (1:d)[prob.estimates <= t]
  result
}

### PART (C) classifier learning with DP SQ algorithm

# ## test
# data.test <- read.csv('hw4testdata.csv')
# sq.centralized(data.test[,names(data.test) != "y"], data.test[,"y"], 0.001, 1)
# sq.local(data.test[,names(data.test) != "y"], data.test[,"y"], 0.5, 1)

data.pums <- read.csv('CaPUMS5Full.csv')
# extract variables data
y.name <- "targetted"
X.names <- names(data.pums)[names(data.pums) != y.name]
X.pums <- data.pums[,X.names]
y.pums <- data.pums[,y.name]

# SQ model results for entire dataset
t.centralized <- 3.558758e-05
t.local <- 0.272

### main results
sqcent <- sq.centralized(X.pums, y.pums, t.centralized, 1)
sqlocal <-sq.local(X.pums, y.pums, t.local, 1)
print(names(X.pums)[sqcent])
print(names(X.pums)[sqlocal])

# calculate false result ratios for a subset S
false.results <- function (X.boot, y.boot, sq, m) {
  n <- length(y.boot)
  # everything is a false positive or true positive if no elements in S
  if (length(sq) == 0) {
    falses <- sum(!as.logical(y.boot)) / n
    return(data.frame(n=n, fp=falses, fn=0, model=m))
  }
  X.boot.sq <- X.boot[,sq, drop=FALSE]
  conj <- apply(X.boot.sq, sum, MARGIN = 1) == length(sq)
  # false positive
  fp.sq <- sum((conj == FALSE) & as.logical(y.boot)) / n
  # false negatives
  fn.sq <- sum((conj == TRUE) & !as.logical(y.boot)) / n
  data.frame(n=n, fp=fp.sq, fn=fn.sq, model=m)
}

# bootstrapping
pums.n <- nrow(X.pums)
results.falses <- data.frame(n=numeric(), fp=numeric(), fn=numeric(), model=character())
# which values of n to look at
nrange <- c(2**(1:log2(100000)), pums.n)
for (i in nrange) {
  boot.rows <- sample(1:pums.n, i, replace=FALSE)
  X.boot <- X.pums[boot.rows,]
  y.boot <- y.pums[boot.rows]
  # centralized
  sqcent <- sq.centralized(X.boot, y.boot, t.centralized, 1)
  results.falses <- rbind(results.falses, false.results(X.boot, y.boot, sqcent, "centralized"))
  # local
  sqloc <- sq.local(X.boot, y.boot, t.local, 1)
  print(i)
  results.falses <- rbind(results.falses, false.results(X.boot, y.boot, sqloc, "local"))
}

library(ggplot2)
# plot false positives and false negatives
plot.falses <- ggplot(results.falses, aes(x=log10(n))) +
  geom_line(aes(y=fp, color=model), linetype="solid") +
  geom_line(aes(y=fn, color=model), linetype="dashed") +
  ylab("false rates") +
  ggtitle("False Positive (solid) & False Negative (dashed) rates on log(n)") +
  theme_bw(); plot.falses
ggsave("falses.jpg", plot.falses)


## bonus 
X.pums.n <- data.frame(!X.pums)
names(X.pums.n) <- paste('n',names(X.pums),sep='')
X.full <- cbind(X.pums, X.pums.n)

sqcentfull <- sq.centralized(X.full, y.pums, t.centralized, 1)
sqlocalfull <- sq.local(X.full, y.pums, t.local, 1)
print(names(X.full)[sqcentfull])
print(names(X.full)[sqlocalfull])