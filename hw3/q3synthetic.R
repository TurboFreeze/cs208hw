## cs208 hw3
## q3.r
##
## Synthetic data release
##
## JH 2019/04/01
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

# Sample with replacement from a vector
bootstrap <- function(x, y, z, n){
  # select indices
  index <- sample(x=1:length(x), size=n, replace=TRUE) 
  list(x=x[index], y=y[index], z=z[index])
}


### Three-way DP histogram
triple.histogram <- function(xdata, ydata, zdata, xlower, xupper, ylower, yupper,
                             zlower, zupper, xnbins=0, ynbins=0, znbins=0, epsilon){
  # x binning
  if(xnbins==0) {
    xlower <- floor(xlower)
    xupper <- ceiling(xupper)
    xbins <- xlower:(xupper+1)    
    xnbins <- length(xbins) -1
    xgranularity <- 1
    xcodebook <- xbins[1:xnbins]
  } else {
    xbins <- seq(from=xlower, to=xupper, length=xnbins+1)
    xgranularity <- (xupper - xlower) / xnbins
    xbins[xnbins+1] <-  xbins[xnbins + 1] + xgranularity
    xcodebook <- xbins[1:xnbins] + 0.5 * xgranularity
  }
  
  # y binning
  if(ynbins==0) {
    ylower <- floor(ylower)
    yupper <- ceiling(yupper)
    ybins <- ylower:(yupper+1)    
    ynbins <- length(ybins) -1
    ygranularity <- 1
    ycodebook <- ybins[1:ynbins]
  } else {
    ybins <- seq(from=ylower, to=yupper, length=ynbins+1)
    ygranularity <- (yupper - ylower) / ynbins
    ybins[ynbins+1] <-  ybins[ynbins + 1] + ygranularity
    ycodebook <- ybins[1:ynbins] + 0.5 * ygranularity
  }
  
  # z binning
  if(znbins==0) {
    zlower <- floor(zlower)
    zupper <- ceiling(zupper)
    zbins <- zlower:(zupper+1)    
    znbins <- length(zbins) -1
    zgranularity <- 1
    zcodebook <- zbins[1:znbins]
  } else {
    zbins <- seq(from=zlower, to=zupper, length=znbins+1)
    zgranularity <- (zupper - zlower) / znbins
    zbins[ynbins+1] <-  zbins[znbins + 1] + zgranularity
    zcodebook <- zbins[1:znbins] + 0.5 * zgranularity
  }
  
  # clip data
  x.clipped <- clip(x=xdata, lower=xlower, upper=xupper)
  y.clipped <- clip(x=ydata, lower=ylower, upper=yupper)
  z.clipped <- clip(x=zdata, lower=zlower, upper=zupper)
  
  # histogram sensitivity and Laplace noise scale
  sensitivity <- 2
  scale <- sensitivity / (epsilon)
  
  sensitive.value <- dp.release <- array(NA, c(xnbins, ynbins, znbins))
  
  # fill in each bin of the three-way histogram
  for(i in 1:xnbins){
    for(j in 1:ynbins){
      for(k in 1:znbins){
        # calculate values that satisfy this bin across all three axis
        sensitive.value[i, j, k] <- sum(x.clipped >= xbins[i] & x.clipped < xbins[i + 1] &
                                          y.clipped >= ybins[j] & y.clipped < ybins[j + 1] &
                                          z.clipped >= zbins[k] & z.clipped < zbins[k + 1])
        # add noise
        dp.release[i, j, k] <- sensitive.value[i, j, k] + rlap(mu=0, s=scale, size=1)
      }
    }
  }
  
  list(release=dp.release, true=sensitive.value, 
       xcodebook=xcodebook, ycodebook=ycodebook, zcodebook=zcodebook)
}


### read data and set parameters
data.pums <- read.csv("MaPUMS5full.csv")
# bootstrap the data
sample.size <- 1000
hist.raw <- bootstrap(x=data.pums$educ, y=data.pums$age, z=data.pums$income, n=sample.size)

# dp histogram configuration
xnbins <- 0
ynbins <- 10
znbins <- 100
xbounds <- c(1, 16)
ybounds <- c(0, 100)
zbounds <- c(0, 100000)

### synthetic data
normalize <- function(x){
  x[x < 0] <- 0
  x <- x / sum(x)
  x
}

# regression results for true values
beta.true <- as.vector(lm(clip(data.pums$income, zbounds[1], zbounds[2])
                          ~ clip(data.pums$educ, xbounds[1], xbounds[2])
                          + clip(data.pums$age, ybounds[1], ybounds[2]))$coef)

# simulate for expectation of beta estimates
nsims <- 100
beta.hats <- matrix(NA, nrow=nsims, ncol=3)
for(i in 1:nsims){
  # make the three-way dp histogram release
  hist.release <- triple.histogram(xdata=hist.raw$x, ydata=hist.raw$y, zdata=hist.raw$z,
                                   xlower=xbounds[1], xupper=xbounds[2], ylower=ybounds[1], yupper=ybounds[2],
                                   zlower=zbounds[1], zupper=zbounds[2],
                                   xnbins=xnbins, ynbins=ynbins, znbins=znbins, epsilon=0.5)
  hist.dp <- hist.release$release
  
  # synthetic data
  syn.prob <- as.vector(normalize(hist.dp))
  syn.xyz <- rmultinom(n=sample.size, prob=syn.prob, size=1)
  syn.educ <- t(syn.xyz) %*% sort(rep(hist.release$xcodebook, dim(hist.dp)[2] * dim(hist.dp)[3]))
  syn.age <- t(syn.xyz) %*% rep(hist.release$ycodebook, dim(hist.dp)[1] * dim(hist.dp)[3])
  syn.income <- t(syn.xyz) %*% sort(rep(hist.release$zcodebook, dim(hist.dp)[1] * dim(hist.dp)[2]))
  
  # linear regression
  output <- lm(syn.income ~ syn.age + syn.educ)
  
  # store calculated regression coefficients
  beta.hats[i,] <- as.vector(output$coef)
}


### Error calculations

calc.mse <- function(actual, estimates) {
  # calculate expectation as mean of beta hats
  estimates.mean <- apply(estimates, 2, mean)
  
  # bias and variance
  bias <- actual - estimates.mean
  var <- apply(apply(estimates, 1, function (x) {t(estimates.mean - x)})^2, 1, mean)
  
  # MSE bias-variance decomposition
  mse <- bias^2 + var
  
  list(biassq=bias^2, var=var, mse=mse)
}

# find dp mse and mse contributions
error.dp <- calc.mse(beta.true, beta.hats)
print(error.dp)

# bootstrap error baseline
beta.bootstrap <- matrix(NA, nrow=nsims, ncol=3)
for(i in 1:nsims){
  # bootstrap the data
  hist.raw <- bootstrap(x=data.pums$educ, y=data.pums$age, z=data.pums$income, n=sample.size)
  # synthetic data
  regression.data <- data.frame(income=clip(hist.raw$z, zbounds[1], zbounds[2]),
                                educ=clip(hist.raw$x, xbounds[1], xbounds[2]),
                                age=clip(hist.raw$y, ybounds[1], ybounds[2]))
  # linear regression and store regression coefficients
  beta.bootstrap[i,] <- as.vector(lm(income ~ educ + age, regression.data)$coef)
}
error.bootstrap <- calc.mse(beta.true, beta.bootstrap)
print(error.bootstrap)

