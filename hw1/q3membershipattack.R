##
## q3membershipattack.r
##
## Membership attack
##
## JH 2019/02/23
##

library(ggplot2)    # using ggplot graphing libraries (optional)

##### Parameters

set.seed(99)              # RNG seed
n <- 100                  # Dataset size
k.trials <- 300      # max number of additional attributes/predicates (m)
n.samples <- 20
n.sims <- 100
delta <- 1 / (10 * n)
# boundary parameters for reconstruction failure
# obtained from reconstruction attack in q2 (see other script)
param.rounding <- 13
param.noisy <- 14
param.subsampling <- 13

##### Data

# read in the data
sample.data <- read.csv('FultonPUMS5sample100.csv')
sample.data.clean <- sample.data[3:17]
attributes.n <- length(names(sample.data.clean))

sample.data.binary <- sample.data.clean[c(1, 5:15)]
binary.names <- names(sample.data.binary)

##### Querying Mechanisms (including defenses)

# standard query of the dataset
query <- function(data, p) {
  # calculate the desired sum
  sum(data) / length(data)
}

# query w/ defense by rounding to nearest multiple of R
query.rounding <- function(data, r) {
  # apply rounding defense
  sum.rounded <- round(sum(data) / r) * r
  sum.rounded / length(data)
}

# query w/ defense by adding Gaussian noise
query.noisy <- function(data, sigma) {
  # apply Gaussian noise as defense
  sum.noisy <- sum(data) + rnorm(1, 0, sigma)
  sum.noisy / length(data)
}

# query w/ defense by subsampling t out of n rows
query.subsampling <- function(data, t) {
  # create subsample T with t out of n rows
  T.rows <- sample(length(data), t)
  # extract data subset from subsample according to specified predicate
  T.data <- data[T.rows]
  # get actual query on this subsample
  sum <- sum(T.data)
  # scale up
  sum.sub <- sum * length(data) / t
  sum.sub / length(data)
}


##### Predicates

# hashing to generate predicates
prime <- 563    # moderately large prime

predicate.single <- function(individual, r.nums) {
  (sum(r.nums * individual) %% prime) %% 2
}

# define predicates for the queries
predicates <- matrix(NA, nrow = k.trials, ncol = n)


# helper function to create data from population
rmvbernoulli <- function(n=1, prob){
  history <- matrix(NA, nrow=n, ncol=length(prob))
  for(i in 1:n){
    x<- rbinom(n=length(prob), size=1, prob=prob)
    x[x==0] <- -1
    history[i,] <- x
  }
  history
}

# test statistic (Dwork et al.)
test.dwork <- function(alice, sample.mean, population.mean){
  sum(alice * sample.mean) - sum(population.mean * sample.mean)
}

# Generate the critical value
generate.critical <- function(null.sims=1000, alpha, fun, population.prob){
  population.mean <- 2 * (population.prob - 0.5)
  hold <- rep(NA,null.sims)
  for(i in 1:null.sims){
    sample <- rmvbernoulli(n=n.samples, prob=population.prob)
    null.alice <- rmvbernoulli(n=1, prob=population.prob)
    sample.mean <- apply(sample, MARGIN=2, FUN=mean)
    hold[i] <- eval(fun(alice=null.alice, sample.mean=sample.mean, population.mean=population.mean))
  }
  null.distribution <- sort(hold, decreasing=TRUE)
  null.distribution[round(alpha * null.sims)]
}


population.data <- read.csv('FultonPUMS5full.csv')
population.probs <- as.vector(apply(population.data[binary.names], 2, mean))
population.means <- 2 * (population.probs - 0.5)

# run simulation against provided query defense mechanism
simulate <- function (query.function, query.param) {  
  predicates <- t(sample.data.binary)
  
  population.probs.all <- population.probs
  population.means.all <- population.means
  
  results <- matrix(NA, nrow = k.trials, ncol = 2)
  
  # loop through number of predicates/attributes
  for (pred.index in 1:k.trials) {
    # generate random numbers for this hash
    r.nums <- sample(0:(prime - 1), attributes.n, replace=TRUE)
    # calculate for particular individual
    pred.temp <- apply(sample.data.clean, MARGIN = 1, predicate.single, r.nums)
    # store the predicate
    predicates <- rbind(predicates, pred.temp)
    population.probs.all <- c(population.probs.all, 0.5)
    population.means.all <- c(population.means.all, 0)
    
    # critical value of the null distribution
    null.critical <- generate.critical(alpha=delta, fun=test.dwork, population.prob=population.probs.all)
    
    history <- vector("numeric", n.sims)
    for(i in 1:n.sims) {
      # Simulate data
      sample <- rmvbernoulli(n=n.samples, prob=population.probs.all)
      sample.mean <- apply(sample, MARGIN=2, FUN=query.function, query.param)
      alice <- sample[1,]
  
      # Conduct hypothesis test
      history[i] <- test.dwork(alice=alice, sample.mean=sample.mean, population.mean=population.means.all)
    }
    
    # fraction identified IN
    accept.frac <- round(100 * mean(history > null.critical)) / 100
    
    # storing results
    results[pred.index, 1] <- nrow(predicates)
    results[pred.index, 2] <- accept.frac
    
    # helps visualize while generating more attributes
    plot(results[, 1], results[, 2], type='l')
  }
  results
}



final.graph <- function(results, defense) {
  # final visualization with ggplot
  ggplot(data.frame(results), aes(x=X1, y=X2)) +
    geom_point(size=0.75) +
    geom_line(alpha=0.5) +
    xlab("Number of queries/attributes") +
    ylab("Fraction of successful membership determination") + 
    ggtitle(paste("Membership Attack Against", defense)) +
    theme_bw()
}



# make simulations and visualizations
simulate.standard <- simulate(query, 0)
final.graph(simulate.standard, "Standard Query w/ No Defense")
dev.copy2pdf(file="./figs/membership.pdf")

simulate.rounding <- simulate(query.rounding, param.rounding)
final.graph(simulate.rounding, "Rounding Defense")
dev.copy2pdf(file="./figs/membershiprounding.pdf")

simulate.noisy <- simulate(query.noisy, param.noisy)
final.graph(simulate.noisy, "Noise Defense")
dev.copy2pdf(file="./figs/membershipnoisy.pdf")

simulate.subsampling <- simulate(query.subsampling, param.subsampling)
final.graph(simulate.subsampling, "Subsampling Defense")
dev.copy2pdf(file="./figs/membershipsubsampling.pdf")

