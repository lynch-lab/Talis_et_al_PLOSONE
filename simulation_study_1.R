
#### Talis et al
# Difficulties in summing distributions for abundance and potential solutions

# simulation study 1

library(latex2exp)

X <- seq(2, 12, 0.1)
sigma <- 0.10
m <- 1000
ensemble_runs <- 10

#########
n <- 10    # redo for all n in c(10,100,1000)
abundance <- array(data=NA, dim=c(ensemble_runs, length(X), n, m))
total.bias.10 <- array(data=NA, dim=c(ensemble_runs, length(X)))

for (k in 1:ensemble_runs){
  for (i in 1:length(X))  # iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance[k,i,j,] <- rlnorm(m,X[i],sigma)  # approximately but not exactly 10% error
    }
  }
  this_ensemble <- abundance[k,,,]
  median.abundance <- apply(this_ensemble, c(1,2), FUN=median) # take the median abundance for each colony
  sum.median <- apply(median.abundance, 1, sum) 
  sum.abundance <- apply(this_ensemble, c(1,3), FUN=sum)
  median.sum <- apply(sum.abundance, 1, median)
  total.bias.10[k,] <- median.sum - sum.median
  # report mean % bias over all mu
  percent.bias.10 <- mean(100*(median.sum-sum.median)/(sum.median))
}
  


###################
n <- 100    # redo for all n in c(10,100,1000)
abundance <- array(data=NA, dim=c(ensemble_runs, length(X), n, m))
total.bias.100 <- array(data=NA, dim=c(ensemble_runs, length(X)))

for (k in 1:ensemble_runs){
  for (i in 1:length(X))  # iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance[k,i,j,] <- rlnorm(m,X[i],sigma)  # approximately but not exactly 10% error
    }
  }
  this_ensemble <- abundance[k,,,]
  median.abundance <- apply(this_ensemble, c(1,2), FUN=median) # take the median abundance for each colony
  sum.median <- apply(median.abundance, 1, sum) 
  sum.abundance <- apply(this_ensemble, c(1,3), FUN=sum)
  median.sum <- apply(sum.abundance, 1, median)
  total.bias.100[k,] <- median.sum - sum.median
  # report mean % bias over all mu
  percent.bias.100 <- mean(100*(median.sum-sum.median)/(sum.median))
}

###################
n <- 1000    # redo for all n in c(10,100,1000)
abundance <- array(data=NA, dim=c(ensemble_runs, length(X), n, m))
total.bias.1000 <- array(data=NA, dim=c(ensemble_runs, length(X)))

for (k in 1:ensemble_runs){
  for (i in 1:length(X))  # iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance[k,i,j,] <- rlnorm(m,X[i],sigma)  # approximately but not exactly 10% error
    }
  }
  this_ensemble <- abundance[k,,,]
  median.abundance <- apply(this_ensemble, c(1,2), FUN=median) # take the median abundance for each colony
  sum.median <- apply(median.abundance, 1, sum) 
  sum.abundance <- apply(this_ensemble, c(1,3), FUN=sum)
  median.sum <- apply(sum.abundance, 1, median)
  total.bias.1000[k,] <- median.sum - sum.median
  # report mean % bias over all mu
  percent.bias.1000 <- mean(100*(median.sum-sum.median)/(sum.median))
  print(k)
}


##################
# plot
library(scales)
plot(log(total.bias.10[1,]),typ="l",ylim=c(-2.5,15),xlab=TeX('$\\mu$'), ylab="log[median(sums) - sum(medians)]", xaxt='n', col=alpha("darkviolet", 0.2)) # 10 colonies
axis(1, at=seq(0, 100, by=20), labels=c(2, 4, 6, 8, 10, 12))
lines(log(total.bias.100[1,]),col=alpha("brown3", 0.2)) # 100 colonies
lines(log(total.bias.1000[1,]),col=alpha("cyan4", 0.2)) # 1000 colonies
legend("topleft", legend=c("n = 1000", "n = 100", "n = 10"), col=c("cyan4", "brown3", "darkorchid1"), lty=1)
for (k in 2:ensemble_runs){
  lines(log(total.bias.10[k,]),col=alpha("darkorchid1", 0.2))
  lines(log(total.bias.100[k,]),col=alpha("brown3", 0.2))
  lines(log(total.bias.1000[k,]),col=alpha("cyan4", 0.2))
}

# average across ensemsbles to get mean line for each n
total.bias.10_mean <- apply(total.bias.10, 2, mean)
total.bias.100_mean <- apply(total.bias.100, 2, mean)
total.bias.1000_mean <- apply(total.bias.1000, 2, mean)
# plolt the mean lines
lines(log(total.bias.10_mean),col="darkorchid1")
lines(log(total.bias.100_mean),col="brown3")
lines(log(total.bias.1000_mean),col="cyan4")

# print mean % bias for each n
#percent.bias.10
#percent.bias.100
#percent.bias.1000

save(total.bias.10, total.bias.100, total.bias.1000, file = "simulation_study_1.Rdata")


