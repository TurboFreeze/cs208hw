## cs208 hw3
## q2.r
##
## DP composition with advanced composition theorem
##
## JH 2019/03/31
##


##### INCLUDING PSILENCE LIBRARY
# Install dependency packages if not current

update_packages <- function(packageList){
  availableRepos <- getCRANmirrors()
  flag <- availableRepos$Country=="USA" & grepl("https",availableRepos$URL,)
  useRepos <- sample(availableRepos$URL[flag],1)
  
  ## install missing packages, and update if newer version available
  for(i in 1:length(packageList)){
    if (!require(packageList[i],character.only = TRUE)){
      install.packages(packageList[i], repos=useRepos)
    }
  }
  
  update.packages(ask = FALSE, dependencies = c('Suggests'), oldPkgs=packageList, repos=useRepos)
}

packagelist <- c("devtools", "jsonlite", "openssl")
update_packages(packagelist)


# Install PSIlence from GitHub
devtools::install_github("privacytoolsproject/PSI-Library", ref="develop") 
#####




# parameter onfiguration
epsilon <- 1
delta <- 10 ^ (-9)
k.range <- 100

# use the PSi library
library("PSIlence")

# basic composition epsilon
calc.basic <- function (epsilon, k) {
  epsilon / k
}

# advanced composition theorem epsilon
calc.advanced <- function (epsilon, k, delta) {
  epsilon / sqrt(2 * k * log(1 / delta))
}

# optimal composition epsilon
calc.optimal <- function (epsilon, k, delta) {
  params <- cbind(rep(1 / k, k), rep(0, k))
  res <- PSIlence:::update_parameters(params=params, hold=0, eps=epsilon, del=delta)
  res[1]
}

# calculate standard deviation for Laplace noise
laplace.std <- function (epsilon0) {
  2 / epsilon0^2
}

# loop through different values of k
stds <- data.frame(k=numeric(), std=numeric(), type=character())
types <- c("basic", "advanced", "optimal")
for (k in 1:k.range) {
  vals <- c(calc.basic(epsilon, k), calc.advanced(epsilon, k, delta),
            calc.optimal(epsilon, k, delta))
  vals <- laplace.std(vals)
  stds.single <- data.frame(k=k, std=vals, type=types)
  stds <- rbind(stds, stds.single)
}

std.basic <- stds[stds$type == "basic",]$std
std.advanced <- stds[stds$type == "advanced",]$std
std.optimal <- stds[stds$type == "optimal",]$std
# find minimum where optimal and advanced both improve on basic compositions
improvement <- min((1:k.range)[std.optimal < std.basic & std.advanced < std.basic])
print(paste("Optimal and Advanced composition improve on basic composition below",
            improvement, "queries"))

# plot
library(ggplot2)
plot.std <- ggplot(stds, aes(x=k)) + 
  geom_line(aes(y=std, color=type), size=1) +
  geom_vline(xintercept=improvement, linetype="dashed") +
  #geom_rect(xmin=0, xmax=improvement, ymin=-10000, ymax=22000, fill="light grey", alpha=0.01) +
  ggtitle("Std. Dev. Against # of Queries for Epsilon Composition") +
  theme_classic()
plot.std

ggsave("std.jpg", plot.std)
