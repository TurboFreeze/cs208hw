##
## q1.r
##
## Exploring Fulton PUMS data to identify variables for reidentification attack 
##
## JH 2019/02/11
##

# read in the data
data <- read.csv('FultonPUMS5.csv')
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
unique(data$puma - data$jpumarow)

# remove variables with single count
results[results$count == 1,] # state and fips have one count
results <- results[results$count != 1,]

# remove extremely private data
results[results$count > 100,]
results <- results[results$count < 100,]
