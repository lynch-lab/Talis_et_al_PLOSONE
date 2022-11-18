#### Talis et al
# Difficulties in summing distributions for abundance and potential solutions

# simulation study performance of alternative solutions

library(truncnorm)

n<-100
X<-seq(2,12,2) 
abundance<-array(data=NA,dim=c(length(X),n,10000))

num_runs <- 100

full_data <- data.frame(correction = NA, X_val = NA, abundance_estimate = NA, order= NA)

# standard log-normal 
total.bias.median.sum <- array(data = NA, dim=c(num_runs, length(X)))

for (k in 1:num_runs) {
  for (i in 1:length(X))  #iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance[i,j,]<-rlnorm(10000,X[i],0.1)  #approximately but not exactly 10% error
    }
  }
  sum.abundance<-apply(abundance, c(1,3), FUN=sum)
  median.sum<-apply(sum.abundance, 1, median) 
  
  total.bias.median.sum[k,] <-100*(median.sum-exp(X)*n)/(exp(X)*n)

  #for (i in 1:length(X)) {
    #uncorrec_data <- data.frame(correction = NA, X_val = NA, abundance_estimate = NA, order= NA)
    #uncorrec_data <- data.frame(correction = "log-normal", X_val = toString(i), abundance_estimate = total.bias.median.sum[k,i], order = 1)
    #full_data <- rbind(full_data, uncorrec_data)
  #}
}

# CORRECTIONS

# correction 1: shifted log-normal
abundance.bias.corrected<-array(data=NA,dim=c(length(X),n,10000))
total.bias.mean.sum.bias.corrected.10 <- array(data = NA, dim=c(num_runs, length(X)))

for (k in 1:num_runs){
  for (i in 1:length(X))  #iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance.bias.corrected[i,j,]<-rlnorm(10000,X[i]-(0.1^2)/2,0.1)  #approximately but not exactly 10% error
    }
  }
  sum.abundance.bias.corrected<-apply(abundance.bias.corrected, c(1,3), FUN=sum)
  mean.sum.bias.corrected<-apply(sum.abundance.bias.corrected, 1, mean) 
  
  total.bias.mean.sum.bias.corrected.10[k,]<-100*(mean.sum.bias.corrected-exp(X)*n)/(exp(X)*n)

  for (i in 1:length(X)) {
    uncorrec_data <- data.frame(correction = NA, X_val = NA, abundance_estimate = NA, order= NA)
    uncorrec_data <- data.frame(correction = "shifted log-normal", X_val = toString(i), abundance_estimate = total.bias.mean.sum.bias.corrected.10[k,i], order = 2)
    full_data <- rbind(full_data, uncorrec_data)
  }
}


####truncated normal ######
abundance.bias.corrected.2<-array(data=NA,dim=c(length(X),n,10000))
total.bias.median.sum.bias.corrected.2.10 <- array(data = NA, dim=c(num_runs, length(X)))

for (k in 1:num_runs){
  for (i in 1:length(X))  #iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance.bias.corrected.2[i,j,]<-rtruncnorm(10000,0,Inf,exp(X[i]),exp(X[i])*0.1)    #approximately but not exactly 10% error
    }
  }
  sum.abundance.bias.corrected.2<-apply(abundance.bias.corrected.2, c(1,3), FUN=sum)
  median.sum.bias.corrected.2<-apply(sum.abundance.bias.corrected.2, 1, mean) 
  
  total.bias.median.sum.bias.corrected.2.10[k,]<-100*(median.sum.bias.corrected.2-exp(X)*n)/(exp(X)*n)

  for (i in 1:length(X)) {
    uncorrec_data <- data.frame(correction = NA, X_val = NA, abundance_estimate = NA, order= NA)
    uncorrec_data <- data.frame(correction = "truncated Normal", X_val = toString(i), abundance_estimate = total.bias.median.sum.bias.corrected.2.10[k,i], order = 3)
    full_data <- rbind(full_data, uncorrec_data)
  }  
}

####censored normal ######
abundance.bias.corrected.3<-array(data=NA,dim=c(length(X),n,10000))
total.bias.median.sum.bias.corrected.3.10 <- array(data = NA, dim=c(num_runs, length(X)))

for (k in 1:num_runs){
  for (i in 1:length(X))  #iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance.bias.corrected.3[i,j,]<-pmax(rnorm(10000,exp(X[i]),exp(X[i])*0.1),0)
    }
  }
  sum.abundance.bias.corrected.3<-apply(abundance.bias.corrected.3, c(1,3), FUN=sum)
  median.sum.bias.corrected.3<-apply(sum.abundance.bias.corrected.3, 1, mean)
  
  total.bias.median.sum.bias.corrected.3.10[k,]<-100*(median.sum.bias.corrected.3-exp(X)*n)/(exp(X)*n)
  
  for (i in 1:length(X)) {
    uncorrec_data <- data.frame(correction = NA, X_val = NA, abundance_estimate = NA, order= NA)
    uncorrec_data <- data.frame(correction = "rectified-before-sum Normal", X_val = toString(i), abundance_estimate = total.bias.median.sum.bias.corrected.3.10[k,i], order = 4)
    full_data <- rbind(full_data, uncorrec_data)
  }
}


#### normal censored post-sum ######
abundance.bias.corrected.4<-array(data=NA,dim=c(length(X),n,10000))
total.bias.median.sum.bias.corrected.4.10 <- array(data = NA, dim=c(num_runs, length(X)))

for (k in 1:num_runs){
  for (i in 1:length(X))  #iterate over mean colony size
  {
    for (j in 1:n) # iterate over the colonies to be summed
    {
      abundance.bias.corrected.4[i,j,]<-rnorm(10000,exp(X[i]),exp(X[i])*0.1)
    }
  }
  sum.abundance.bias.corrected.4<-apply(abundance.bias.corrected.4, c(1,3), FUN=sum)
  sum.abundance.bias.corrected.4<-pmax(sum.abundance.bias.corrected.4,0)
  median.sum.bias.corrected.4<-apply(sum.abundance.bias.corrected.4, 1, mean) 
  
  total.bias.median.sum.bias.corrected.4.10[k,]<-100*(median.sum.bias.corrected.4-exp(X)*n)/(exp(X)*n)
  
  for (i in 1:length(X)) {
    uncorrec_data <- data.frame(correction = NA, X_val = NA, abundance_estimate = NA, order= NA)
    uncorrec_data <- data.frame(correction = "rectified-after-sum Normal", X_val = toString(i), abundance_estimate = total.bias.median.sum.bias.corrected.4.10[k,i], order = 4)
    full_data <- rbind(full_data, uncorrec_data)
  }
}

#########Now compare everything ##########
#plot(total.bias.sum.medians,typ="l",lwd=2,ylim=c(-800,800), xlab="mean(logged-abundance)", ylab="total abundance - exp(X)*n", xaxt='n')
#axis(1, at=seq(0, 100, by=20), labels=c(2, 4, 6, 8, 10, 12))
#lines(total.bias.median.sum,col="red",lwd=2)
#lines(total.bias.sum.median.bias.corrected.3.10,col="yellow",lwd=2)
#lines(total.bias.median.sum.bias.corrected.3.10,col="pink",lwd=2)
#lines(total.bias.sum.median.bias.corrected.2.10,col="blue",lwd=2)
#lines(total.bias.median.sum.bias.corrected.2.10,col="purple",lwd=2)
#lines(total.bias.sum.mean.bias.corrected.10,col="orange",lwd=2)
#lines(total.bias.mean.sum.bias.corrected.10,col="green",lwd=2)

# library
library(ggplot2)
library(forcats)
library(dplyr)

full_data <- full_data[-c(1),]

save(full_data, file = "fig_2_mean_ensemble.Rdata")

full_data <- mutate(full_data, correction = fct_reorder(correction, order))

full_data$correction <- factor(full_data$correction, levels = c('shifted log-normal', 'truncated Normal', 'rectified-before-sum Normal', 'rectified-after-sum Normal'))

library(latex2exp)

load("fig_2_mean_ensemble.Rdata")

# grouped boxplot
ggplot(full_data, aes(x=X_val, y=abundance_estimate, fill=factor(correction))) + 
  geom_boxplot() + theme_minimal() +
  #coord_cartesian(ylim = c(-0.6, 0.6)) +
  scale_x_discrete(labels = c(2, 4, 6, 8, 10, 12)) +
  xlab(TeX('$\\mu$')) +
  ylab("% bias") +
  #ylab("(total abundance - exp(X)*n)/(exp(X)*n)*100") +
  guides(fill=guide_legend(title="Correction"))



print(c(mean(total.bias.median.sum), sd(total.bias.median.sum)))
print(c(mean(total.bias.mean.sum.bias.corrected.10), sd(total.bias.mean.sum.bias.corrected.10)))
print(c(mean(total.bias.median.sum.bias.corrected.2.10), sd(total.bias.median.sum.bias.corrected.2.10)))
print(c(mean(total.bias.median.sum.bias.corrected.3.10), sd(total.bias.median.sum.bias.corrected.3.10)))
print(c(mean(total.bias.median.sum.bias.corrected.4.10), sd(total.bias.median.sum.bias.corrected.3.10)))
