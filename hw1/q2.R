##
## q2.r
##
## Regression based attack against various defense mechanisms
##
## JH 2019/02/21
##

library(ggplot2)    # using ggplot graphing libraries (optional)

##### Parameters

set.seed(123)       # RNG seed
n <- 100            # Dataset size
k.trials <- 2 * n   # Number of queries
exp.n <- 10          # Number of experiments at each step

##### Data

# read in the data
sample.data <- read.csv('FultonPUMS5sample100.csv')

# get sensitive data (USCITIZEN)
sensitive.var <- "uscitizen"
sensitive.data <- sample.data[, sensitive.var]


##### Querying Mechanisms (including defenses)

# standard query of the dataset
query <- function(data, pred) {
  # extract data subset according to specified predicate
  subset <- data[pred]
  # calculate desired sum
  sum <- sum(subset)
  sum
}

# query w/ defense by rounding to nearest multiple of R
query.rounding <- function(data, pred, r) {
  # get actual query
  sum <- query(data, pred)
  # apply rounding defense
  sum.rounded <- round(sum / r) * r
  sum.rounded
}

# query w/ defense by adding Gaussian noise
query.noisy <- function(data, pred, sigma) {
  # get actual query
  sum <- query(data, pred)
  # apply Gaussian noise as defense
  sum.noisy <- sum + rnorm(1, 0, sigma)
  sum.noisy
}

# query w/ defense by subsampling t out of n rows
query.subsampling <- function(data, pred, t) {
  # create subsample T with t out of n rows
  T.rows <- sample(length(data), t)
  # extract data subset from subsample according to specified predicate
  T.data <- data[T.rows]
  # get actual query on this subsample
  sum <- query(T.data, pred[T.rows])
  # scale up
  sum.sub <- sum * length(data) / t
  sum.sub
}


##### Define Predicates

# define predicates for the 2n queries
predicates <- matrix(NA, nrow = k.trials, ncol = n)

# helper function to determine membership in pred matrix
pred.membership <- function(row.temp) {
  bool.membership <- FALSE
  for (i in 1:200) {
    bool.membership <- bool.membership ||
      isTRUE(all.equal(row.temp, predicates[i,]))
  }
  bool.membership
}

# generate 2n predicates
for (pred.index in 1:k.trials) {
  # generate possible predicate (boolean array)
  pred.temp <- sample(c(TRUE, FALSE), n, replace = TRUE)
  # keep generating until unique found
  # (inefficient but workable for this size)
  while (pred.membership(pred.temp)) {
    pred.temp <- sample(c(TRUE, FALSE), n, replace = TRUE)
  }
  # store the predicate
  predicates[pred.index,] <- pred.temp
}


#### Main Reconstruction Attack function
# takes a defense query mechanism and name of defense
# makes the attack and plots the results
reconstruction <- function (query.function, defense.name) {
  results <- data.frame(param=integer(), rmse=numeric(), frac=numeric())
  # loop through different parameter values
  for (param.value in 1:n) {
    # make 10 experiments each time, storing each RMSE and fraction of success
    experiments.rmse <- vector("numeric", exp.n)
    experiments.frac <- vector("numeric", exp.n)
    for (exp.index in 1:exp.n) {
      squared.errors <- c()
      history <- matrix(NA, nrow = k.trials, ncol = (n + 1))
      
      ##### Run Queries
      # make the 2n queries using the generated predicates
      for (pred.index in 1:k.trials) {
        # extract predicate
        pred.current <- predicates[pred.index,]
        # get standard query and query with defense
        q.standard <- query(sensitive.data, pred.current)
        q.defense <- query.function(sensitive.data, pred.current, param.value)
        # save this query (result w/ defense + predicate) to history table
        history[pred.index,] <- c(q.defense, as.numeric(pred.current))
        # calculate and store squared errors
        squared.errors <- c(squared.errors, (q.defense - q.standard) ** 2)
      }
      # calculate RMSE for this experiment
      experiments.rmse[exp.index] <- sqrt(sum(squared.errors))
      
      # convert into dataframe
      xnames <- paste("x", 1:n, sep="")
      varnames <- c("y", xnames)
      release.data <- as.data.frame(history)
      names(release.data) <- varnames
      
      ##### Regression Attack
      # attack formula
      attack.formula <- paste(xnames, collapse=" + ")
      attack.formula <- paste("y ~ ", attack.formula, "-1")
      attack.formula <- as.formula(attack.formula)
      attack.output <- lm(attack.formula, data=release.data)
      # regression coefficient estimates
      attack.estimates <- attack.output$coef
      # estimates rounded to give reconstruction predictions
      attack.reconstructed <- round(attack.estimates)
      # calculate fraction successfully reconstructed
      experiments.frac[exp.index] <-
        sum(attack.reconstructed == sensitive.data) / n
    }
    # append the avg of the 10 experiments for this particular parameter value
    results <- rbind(results,
                     c(param.value, mean(experiments.rmse), mean(experiments.frac)))
  }
  
  names(results) <- c("param", "rmse", "frac")
  
  
  ##### Visualization of results
  # plot the result
  ggplot(results, aes(x=rmse, y=frac)) + 
    geom_point(aes(color=param)) + 
    ylab("Fraction of Successful Reconstruction") +
    xlab("RMSE") +
    ggtitle(paste("Reconstruction Attack against", defense.name)) + 
    theme_bw()
  # plot(results$frac, results$rmse) # use this instead if no ggplot
}

# make the calls to the reconstruction attack function
reconstruction(query.rounding, "Rounding Defense")
reconstruction(query.noisy, "Noise Defense")
reconstruction(query.subsampling, "Subsampling Defense")
