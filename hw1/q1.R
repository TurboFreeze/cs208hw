##
## q1.r
##
## Exploring Fulton PUMS data to identify variables for reidentification attack 
##
## JH 2019/02/11
##

# read in the data
data <- read.csv('FultonPUMS5full.csv')
# summary of data
head(data)
nrow(data)
# variables
print(paste(names(data), collapse=", "))
print(length(names(data)))

# see how many unique values of each variable
vars <- c()
counts <- c()
for (n in names(data)) {
  vars <- c(vars, n)
  counts <- c(counts, nrow(unique(data[n])))
}
results <- data.frame(var=vars, count=counts)

# look at location identifiers
table(data$puma)
# unique(data$puma - data$jpumarow)
# unique(data$puma - data$X)
# # remove duplicates
# results <- results[!results$var %in% c("X", "jpumarow"),]

# remove variables with single count
results[results$count == 1,] # state and fips have one count
results <- results[results$count != 1,]


# plot distributions
jpeg('./figs/exploredist.jpg', width=1000, height=300)
par(mfrow=c(1,3))
barplot(table(data$income), xaxt='n', ylab='', main="income Distribution")
axis(side=1)
barplot(table(data$age), ylab='', main="age Distribution")
barplot(table(data$educ), ylab='', main="educ Distribution")
par(mfrow=c(1,1))
dev.off()


# collision probabilities
collision <- function(vardata) {
  sum((table(vardata) / nrow(data))^2)
}

p.collisions <- c()
for (v in as.character(results$var)) {
  p.collisions <- rbind(p.collisions, c(v,collision(data[v])))
}
p.collisions

# find average size of PUMA
avg.puma <- round(20 * nrow(data) / length(unique(data$puma)))

# calculate joint probability of collision
collision.joint <- function (vars) {
  joint <- 1
  for (v in vars) {
    joint <- joint * collision(data[v])
  }
  joint
}

# calculate percentage of average PUMA uniquely identifiable
identifiable <- function (p.collision) {
  (1 - p.collision)^avg.puma
}


# calculate results with proposed list of identifying variables
varlist <- c("income", "age", "educ")
identifiable(collision.joint(varlist))
varlist <- c(varlist, "sex", "married")
identifiable(collision.joint(varlist))
varlist <- c(varlist, "black")
identifiable(collision.joint(varlist))
varlist <- c(varlist, "employed")
identifiable(collision.joint(varlist))
