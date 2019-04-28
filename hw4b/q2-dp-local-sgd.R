## cs208 hw4b
## q2-dp-local-sgd.r
##
## local differentially private stochastic gradient descent (DP-SGD)
##
## JH 2019/04/26
##

library(ggplot2)
library(scales)

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

# negative log-likelihood of logistic regression
calc.logreglik <- function(b, data){ 
  y <- as.integer(data[1])
  x <- as.integer(data[2])
  
  # systematic component
  pi <- 1 / (1 + exp(-b[1] - b[2]*x))
  
  # stochastic component
  llik <- y * log(pi) + (1 - y) * log(1 - pi)
  -llik
}

# Differentially Private Gaussian release
mechanism.gaussian <- function(size=1, sensitivity, epsilon, delta){
  scale <- sensitivity * sqrt(2 * log(1.25 / delta)) / epsilon
  noise <- rnorm(n=size, mean=0, sd=scale)
  noise
}

# helper function for calculating RMSE
rmse <- function (actual, release) {
  sqrt(sum((actual - release) ^ 2))
}

# logistic regression prediction
predict.logreg <- function (data.x, coefs) {
  round(1 / (1 + exp(-(coefs[2] * data.x + coefs[1]))))
}

# calculate the gradient
calc.grad <- function(B, C, theta, fun=calc.logreglik) {
  dx <- 0.0001
  out1 <-	eval(fun(b=theta, data=B))
  out2 <- eval(fun(b=theta + c(0,dx), data=B))
  out3 <- eval(fun(b=theta + c(dx,0), data=B))
  
  Del.1 <- (out3 - out1)/dx
  Del.1 <- clip(Del.1,-C,C)
  mean.Del.1 <- mean(Del.1)
  
  
  Del.2 <- (out2 - out1)/dx
  Del.2 <- clip(Del.2,-C,C)
  mean.Del.2 <- mean(Del.2)
  
  c(mean.Del.1,mean.Del.2)
}

# calculate local dp gradient
calc.grad.dp.local <- function (B, C, theta, sensitivity, epsilon, delta) {
  # clipped true gradient
  grad <- calc.grad(B, C, theta)
  # noisy gradient
  grad.noisy <- grad + mechanism.gaussian(2, sensitivity, epsilon, delta)
  grad.noisy
}

# read the data
pums.data <- read.csv("MaPUMS5full.csv")
pums.data.subset <- pums.data[c("married", "educ")]

n <- nrow(pums.data) # dataset size

# parameters
delta <- 10 ** -6                             # privacy delta
epsilons <- 10 ** (seq(-2, 1, 0.2))          # privacy epsilons to try
steps <- round(sqrt(n))                       # number of steps
nu <- 0.01                                    # step size
clipper <- 10                                 # value to clip to
batch.size <- round(sqrt(n))                  # batch size
ind.limit <- ceiling(steps * batch.size / n)  # limit for individual to be in
theta.start <- c(0, 0)                        # Starting estimate

# locally differentially private stochastic gradient descent (local dp-sgd)
dp.sgd.local <- function (data.shuffle, epsilon) {
  theta <- theta.start
  # loop through gradient descent steps
  for (i in 1:steps) {
    # extract batch 
    batch.start <- ((i - 1) * batch.size + 1)
    if (i < batch.size){
      batch.stop <- i * batch.size
    } else {
      batch.stop <- n
    }
    batch.data <- data.shuffle[batch.start:batch.stop, ]
    
    #Del <- calcgradient(B, C, theta, fun=calcllik) + gaussianReleaseNoise(size=2, sensitivity=2*C/L, epsilon=epsilon*L/2, delta=1e-5)
    
    # local differentially private estimates of gradients
    Del.locals <- apply(batch.data, MARGIN=1, calc.grad.dp.local, clipper, theta, 2 * clipper, epsilon * steps / 2, delta)

    # server averaging of local estimates
    Del <- apply(Del.locals, mean, MARGIN=1)
    
    # update coefficient estimate
    theta <- theta - (Del * nu)
  }
  theta
}


# shuffle the data
pums.data.shuffle <- pums.data.subset[sample(1:nrow(pums.data.subset)),]

epsilon.n <- length(epsilons)
rmses <- vector("numeric", epsilon.n)
classerror <- vector("numeric", epsilon.n)

# actual logistic regression results
coefs.actual <- coef(glm(married ~ educ, family="binomial", data=pums.data.subset))

# calculate rmse of local dp sgd with different epsilon
for (i in 1:epsilon.n) {
  # differentially private local sgd coefficient estimates
  coefs.dp <- dp.sgd.local(pums.data.shuffle, epsilons[i])
  # RMSE
  rmses[i] <- rmse(coefs.actual, coefs.dp)
  # classification error
  predictions.dp <- predict.logreg(pums.data.shuffle$educ, coefs.dp)
  classerror[i] <- sum(abs(predictions.dp - pums.data.shuffle$married)) / n
}

results <- data.frame(epsilon=epsilons, rmse=rmses, cerror=classerror)
# plotting RMSEs
plot.rmse <- ggplot(results, aes(x=epsilon)) +
  geom_line(aes(y=rmse), color="red") + # rmse
  scale_x_continuous(trans='log10', labels=number_format(accuracy = 0.001), breaks = trans_breaks('log10', function(x) 10^x)) +
  ggtitle("RMSE of local DP-SGD for varying epsilon") +
  theme_bw(); plot.rmse
ggsave("rmse1.jpg", plot.rmse)

plot.cerror <- ggplot(results, aes(x=epsilon)) +
  geom_line(aes(y=cerror), color="blue") + # classification error
  scale_x_continuous(trans='log10', labels=number_format(accuracy = 0.001), breaks = trans_breaks('log10', function(x) 10^x)) +
  ggtitle("Classification Error of local DP-SGD for varying epsilon") +
  theme_bw(); plot.cerror
ggsave("cerror1.jpg", plot.cerror)


