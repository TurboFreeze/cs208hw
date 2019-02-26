##
## q2reconstructionattack.r
##
## Regression-based reconstruction attack against various defense mechanisms
##
## JH 2019/02/21
##

library(ggplot2)    # using ggplot graphing libraries (optional)

##### Parameters

set.seed(99)          # RNG seed
n <- 100              # Dataset size
k.trials <- 2 * n     # Number of queries
exp.n <- 10           # Number of experiments at each step

##### Data

# read in the data
sample.data <- read.csv('FultonPUMS5sample100.csv')
sample.data.clean <- sample.data[3:17]
attributes.n <- length(names(sample.data.clean))

# get sensitive data (USCITIZEN)
sensitive.var <- "uscitizen"
sensitive.data <- sample.data[, sensitive.var]


##### Querying Mechanisms (including defenses)

# standard query of the dataset
query <- function(data, pred, p) {
  # extract data subset according to specified predicate
  subset <- data[as.logical(pred)]
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
  sum <- query(data[T.rows], pred[T.rows])
  # scale up
  sum.sub <- sum * length(data) / t
  sum.sub
}


##### Define Predicates

# hashing to generate predicates
prime <- 563    # moderately large prime

predicate.single <- function(r.nums, individual) {
  (sum(r.nums * individual) %% prime) %% 2
}

# define predicates for the 2n queries
predicates <- matrix(NA, nrow = k.trials, ncol = n)

# generate 2n predicates
for (pred.index in 1:k.trials) {
  # generate random numbers for this hash
  r.nums <- sample(0:(prime - 1), attributes.n, replace=TRUE)
  # calculate for particular individual
  pred.temp <- apply(sample.data.clean, MARGIN = 1, predicate.single, r.nums)
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
      attack.reconstructed <- as.numeric(attack.estimates > 0.5) #round(attack.estimates)
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
  # type to distinguish between successful (1) or not (0)
  results$success <- as.numeric(results$frac > 0.5)
  # plot the result
  ggplot(results, aes(x=rmse, y=frac)) + # scatter plot
    geom_point(aes(color=param)) + 
    # trend line
    geom_line(aes(color=param), alpha=0.3) +
    # success threshold
    geom_line(aes(y=0.5), alpha=0.5, size=0.5, linetype=2) +
    # labels and title
    ylab("Fraction of Successful Reconstruction") +
    xlab("RMSE") +
    ggtitle(paste("Reconstruction Attack against", defense.name)) + 
    # legend formatting
    scale_color_continuous(name="Parameter\nValue", low="#56B1F7", high="#132B43") +
    theme_bw()
  #plot(results$frac, results$rmse) # use this instead if no ggplot
}

# make the calls to the reconstruction attack function
reconstruction(query.rounding, "Rounding Defense")
dev.copy2pdf(file="./figs/attackrounding.pdf")
reconstruction(query.noisy, "Noise Defense")
dev.copy2pdf(file="./figs/attacknoise.pdf")
reconstruction(query.subsampling, "Subsampling Defense")
dev.copy2pdf(file="./figs/attacksubsampling.pdf")
