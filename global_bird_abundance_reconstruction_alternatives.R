
#### Talis et al
# Difficulties in summing distributions for abundance and potential solutions


### Global bird abundance simulations
# (Callaghan et al data reconstruction & performance of alternative solutions)

library(ggplot2)
library(gridExtra)
library(rmarkdown)
library(fitdistrplus)
library(gtable)
library(draw)
library(truncnorm)
library(moments)



# load callaghan data
bird_abundances <- read.csv("pnas.2023170118.sd01.csv")

colnames(bird_abundances) <- c("CommonName", "ScientificName", "Order", "Family", "LowerCI", "AbundanceEstimate", "UpperCI", "RangeAdjusted", "TrainingSpecies")

paged_table(bird_abundances)

bird_abundances$sdnorm <- NA
bird_abundances$sdlog <- NA
bird_abundances$meanlog <- NA


### Fitting log-normal distributions to given confidence intervals
LL<-691718494	
median<-1618744682	
UL<-3795135209
log10.LL<-log10(LL)
log10.median<-log10(median)
log10.UL<-log10(UL)
est.sd<-mean(c(log10.median-log10.LL,log10.UL-log10.median))/1.96
simmed.vals<-10^(rnorm(1000,mean=log10.median,sd=est.sd))
fitted.params<-fitdistr(simmed.vals,"lognormal")
hist(simmed.vals,breaks=30)
hist(rlnorm(1000,fitted.params$estimate[1],fitted.params$estimate[2]),add=T,col=rgb(1,0,0,0.5),breaks=30)


rtruncated_norm <- function(n, mu, sigma, low, high) {
  # find quantiles that correspond the the given low and high levels.
  p_low <- pnorm(low, mu, sigma)
  p_high <- pnorm(high, mu, sigma)
  
  # draw quantiles uniformly between the limits and pass these
  # to the relevant quantile function.
  qnorm(runif(n, p_low, p_high), mu, sigma)
}



num_species <- length(bird_abundances$CommonName)
num_samples <- 10000

bird_normals <- array(data=NA, dim=c(num_species, num_samples))
bird_abundances$Our_abundance_estimate <- NA

bird_lognormals_raw <- array(data=NA, dim=c(num_species, num_samples))
bird_lognormals_sort <- array(data=NA, dim=c(num_species, num_samples))
bird_lognormals_shift <- array(data=NA, dim=c(num_species, num_samples))
bird_lognormals_trunc <- array(data=NA, dim=c(num_species, num_samples))
bird_lognormals_recbefore <- array(data=NA, dim=c(num_species, num_samples))
bird_lognormals_recafter <- array(data=NA, dim=c(num_species, num_samples))


sd_norm <- function(mu, sd){
  sd_n <- sqrt((exp(sd^2)-1)*(exp(2*mu + sd^2)))
  return (sd_n)
}

# fit a normal to log base 10 (confidence intervals)
for (species in 1:num_species){
  
  #species <- 1
  
  median_abundance <- bird_abundances[species,]$AbundanceEstimate
  lowerCI_abundance <- bird_abundances[species,]$LowerCI
  upperCI_abundance <- bird_abundances[species,]$UpperCI
  
  log10_median <- log10(median_abundance)
  if (log10_median == -Inf){
    log10_median <- 0
  } else if (log10_median == Inf){
    log10_median <- 10
  }
  log10_lowerCI <- log10(lowerCI_abundance)
  if (log10_lowerCI == -Inf){
    log10_lowerCI <- 0
  }
  log10_upperCI <- log10(upperCI_abundance)
  if (log10_upperCI == Inf){
    log10_lowerCI <- 10
  }
  
  bird_abundances[species,]$sdnorm <- mean(c(log10_median-log10_lowerCI,log10_upperCI-log10_median))/1.96
  bird_normals[species,] <- 10^(rnorm(num_samples, mean=log10_median, sd=bird_abundances[species,]$sdnorm))
  
  lognormal_params <- fitdistr(bird_normals[species,], "lognormal")
  bird_abundances[species,]$meanlog <- as.numeric(lognormal_params$estimate[1])
  bird_abundances[species,]$sdlog <- as.numeric(lognormal_params$estimate[2])
  
  # raw
  bird_lognormals_raw[species,] <- rlnorm(num_samples, bird_abundances[species,]$meanlog, bird_abundances[species,]$sdlog)
  bird_abundances[species,]$Our_abundance_estimate <- median(bird_lognormals_raw[species,])
  
  # sorted
  bird_lognormals_sort[species,] <- sort(bird_lognormals_raw[species,])
  
  # shifted
  bird_lognormals_shift[species,] <- rlnorm(num_samples, bird_abundances[species,]$meanlog - (bird_abundances[species,]$sdlog^2)/2, bird_abundances[species,]$sdlog)
  
  sd_n <- sd_norm(bird_abundances[species,]$meanlog, bird_abundances[species,]$sdlog)
  
  # truncated normal
  bird_lognormals_trunc[species,] <- rtruncated_norm(num_samples, exp(bird_abundances[species,]$meanlog), sd_n, 0, Inf)
  
  # rectified (before summing) normal
  bird_lognormals_recbefore[species,] <- pmax(0, rnorm(num_samples, exp(bird_abundances[species,]$meanlog), sd_n))
  
  # normal (to be rectified after summing)
  bird_lognormals_recafter[species,] <- rnorm(num_samples, exp(bird_abundances[species,]$meanlog), sd_n)
}


df_ex <- data.frame(draw = rlnorm(10000, bird_abundances[1,]$meanlog, bird_abundances[1,]$sdlog)) 

plot3 <- ggplot(df_ex, aes(draw)) + geom_histogram(color="black", fill="darkgray", bins = 100) + geom_density(color = "black", fill = "blue") + xlab("Number of individual birds (linear)")
plot3


hist(log(rlnorm(10000, log(bird_abundances[1,]$AbundanceEstimate), bird_abundances[1,]$sdlog)), breaks = 100)



hist(log(bird_abundances$Our_abundance_estimate, base = 10), breaks = 100)



for (row in which(bird_abundances$Our_abundance_estimate == 0)){ # add a constant 1 for those species predicted to have 0 abundance
  bird_abundances[row,]$Our_abundance_estimate <- 1
}
p2A <- ggplot(bird_abundances, aes(Our_abundance_estimate)) + geom_histogram(binwidth = .20, color="black", fill="darkgray") + ylim(0, 800) + xlab("Species Abundance Distribution") + ylab("Number of species") + coord_flip() + scale_x_log10(limits = c(0.1, 1e10), breaks = c(100, 1e6, 1e9), label = c("100","100,000","100,000,000"))
p2A

ggsave("callaghanFig2A.eps")


# Fig 2B part 1
ind <- which(bird_abundances$CommonName == "Ring-billed Gull")
df_2B_1 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_1 <- ggplot(df_2B_1, aes(draw)) + geom_density(bins = 100, color="black", fill="dimgray") + xlab("") + ylab("") + ylim(0, 0.6) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Green Heron")
df_2B_2 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_2 <- ggplot(df_2B_2, aes(draw)) + geom_density(bins = 100, color="black", fill="lightgray") + xlab("") + ylab("") + ylim(0, 0.6) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Northern Wheatear")
df_2B_3 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_3 <- ggplot(df_2B_3, aes(draw)) + geom_density(bins = 100, color="black", fill="dimgray") + xlab("") + ylab("") + ylim(0, 0.7) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Ashy Prinia")
df_2B_4 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_4 <- ggplot(df_2B_4, aes(draw)) + geom_density(bins = 100, color="black", fill="lightgray") + xlab("") + ylab("") + ylim(0, 0.6) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Osprey")
df_2B_5 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_5 <- ggplot(df_2B_5, aes(draw)) + geom_density(bins = 100, color="black", fill="dimgray") + xlab("") + ylab("") + ylim(0, 1.2) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Acorn Woodpecker")
df_2B_6 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_6 <- ggplot(df_2B_6, aes(draw)) + geom_density(bins = 100, color="black", fill="lightgray") + xlab("") + ylab("") + ylim(0, 0.6) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Yellow-tailed Black-Cockatoo")
df_2B_7 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_7 <- ggplot(df_2B_7, aes(draw)) + geom_density(bins = 100, color="black", fill="dimgray") + xlab("") + ylab("") + ylim(0, 0.6) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))

ind <- which(bird_abundances$CommonName == "Midget Flowerpecker")
df_2B_8 <- data.frame(draw = bird_lognormals_raw[ind,])
p2B_8 <- ggplot(df_2B_8, aes(draw)) + geom_density(bins = 100, color="black", fill="lightgray") + xlab("Number of individual birds (millions)") + ylab("") + ylim(0, 0.6) + scale_x_log10(limits = c(10, 1e13), breaks = c(1e3, 1e6, 1e9, 1e12), label = c("0","1","1,000","1,000,000"))


grid1 <- grid.arrange(p2B_1, p2B_2, p2B_3, p2B_4, nrow = 4, ncol = 1)
grid1

ggsave("callaghanFig2B_1.eps", plot = grid1)


# Fig 2B part 2
grid2 <- grid.arrange(p2B_5, p2B_6, p2B_7, p2B_8, nrow = 4, ncol = 1)
grid2

ggsave("callaghanFig2B_2.eps", plot = grid2)


median(bird_abundances$Our_abundance_estimate) # (their estimate: 450,000)
mean(bird_abundances$Our_abundance_estimate)/1000000 # in millions (their estimate: 5.2)
sum(bird_abundances$Our_abundance_estimate)/1000000000 # sum of medians in billions:



#### Raw (Uncorrected)
total_dist_raw <- c()

for (draw in 1:num_samples){
  
  total_dist_raw[draw] <- sum(bird_lognormals_raw[,draw])
  
}


df_2C_raw <- data.frame(sum = total_dist_raw)
sd_raw <- round(sd(total_dist_raw)/1000000000,3) # in billions

p2C_raw <- ggplot(df_2C_raw, aes(sum)) + geom_density(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") + geom_vline(xintercept = 392000000000, color = "black") + geom_vline(xintercept = 428000000000, color = "blue") # + geom_label(aes(x=392000000000, label="392 billion", y=4), colour="black") + geom_label(aes(x=428000000000, label="428 billion", y=3.5), colour="blue") + geom_label(aes(x=5e12, label=paste("sd =",sd_raw,"billion"), y=2.5), colour="black")
p2C_raw

ggsave("raw_lognormal.eps")

l2C_raw <- ggplot(df_2C_raw, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_discrete(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
l2C_raw

ggsave("raw_lognormal_linear.eps")

median(total_dist_raw)/1000000000 # in billions
mean(total_dist_raw)/1000000000 # in billions


#### Sorted
total_dist_sort <- c()

for (draw in 1:num_samples){
  
  total_dist_sort[draw] <- sum(bird_lognormals_sort[,draw])
  
}


df_2C_sort <- data.frame(sum = total_dist_sort)
sd_sort <- round(sd(total_dist_sort)/1e12,3) # in trillions

p2C_sort <- ggplot(df_2C_sort, aes(sum)) + geom_density(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") + geom_vline(xintercept = 50000000000, color = "black") + geom_vline(xintercept = 428000000000, color = "blue") # + geom_label(aes(x=50000000000, label="50 billion", y=0.60), colour="black") + geom_label(aes(x=428000000000, label="428 billion", y=0.60), colour="blue") + geom_label(aes(x=5e12, label=paste("sd =",sd_sort,"trillion"), y=0.35), colour="black")
p2C_sort

ggsave("sorted_lognormal.eps")



l2C_sort <- ggplot(df_2C_sort, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_discrete(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
l2C_sort

l2C_sort <- ggplot() + 
  geom_step(aes(df_2C_sort$sum, y =..density..),
            stat = 'bin',breaks = dhist(df_2C_sort$sum, nbins = 30),
            position = "dodge", color = 'red')

breaks = c(0, 0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13)

l2C_sort <- ggplot(df_2C_sort, aes(sum)) +
  geom_histogram(aes(y=..density..),
                 color="black", fill="grey40", breaks=breaks) +
  scale_x_continuous(breaks=breaks, limits = c(0, 1e13))
l2C_sort

ggsave("sorted_lognormal_linear.eps")



median(total_dist_sort)/1000000000 # in billions (their estimate: 50)
mean(total_dist_sort)/1000000000 # in billions (their estimate: 428)

(max(total_dist_sort) - min(total_dist_sort))/1000000000



#### Shifted Log-normal
total_dist_shift <- c()

for (draw in 1:num_samples){
  
  total_dist_shift[draw] <- sum(bird_lognormals_shift[,draw])
  
}



df_2C_shift <- data.frame(sum = total_dist_shift)
sd_shift <- round(sd(total_dist_shift)/1000000000,3) # in billions

p2C_shift <- ggplot(df_2C_shift, aes(sum)) + geom_density(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") + geom_vline(xintercept = 47000000000, color = "black") + geom_vline(xintercept = 50000000000, color = "blue")  # + geom_label(aes(x=47000000000, label="47 billion", y=6.1), colour="black") + geom_label(aes(x=50000000000, label="50 billion", y=5.5), colour="blue") + geom_label(aes(x=5e12, label=paste("sd =",sd_shift,"billion"), y=3.25), colour="black")
p2C_shift

ggsave("shifted_lognormal.eps")

l2C_shifted <- ggplot(df_2C_shift, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_discrete(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
l2C_shifted

ggsave("shifted_lognormal_linear.eps")


median(total_dist_shift)/1000000000 # in billions
mean(total_dist_shift)/1000000000 # in billions



#### Truncated Normal
total_dist_trunc <- c()

for (draw in 1:num_samples){
  
  total_dist_trunc[draw] <- sum(bird_lognormals_trunc[,draw])
  
}



df_2C_trunc <- data.frame(sum = total_dist_trunc)
sd_trunc <- round(sd(total_dist_trunc)/1000000000,3) # in billions

p2C_trunc <- ggplot(df_2C_trunc, aes(sum)) + geom_density(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 94300000000, color = "black") + geom_vline(xintercept = 94400000000, color = "blue") # + geom_label(aes(x=94300000000, label="94.3 billion", y=30), colour="black") + geom_label(aes(x=94400000000, label="94.4 billion", y=27), colour="blue") + geom_label(aes(x=5e12, label=paste("sd =",sd_trunc,"billion"), y=20), colour="black")
p2C_trunc

ggsave("truncated_normal.eps")


l2C_trunc <- ggplot(df_2C_trunc, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_discrete(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
l2C_trunc

ggsave("truncated_normal_linear.eps")


median(total_dist_trunc)/1000000000 # in billions
mean(total_dist_trunc)/1000000000 # in billions

means <- c()
medians <- c()
diffs <- c()
diffs_sym <- c()
for (species in 1:num_species){
  means[species] <- mean(bird_lognormals_trunc[species,])
  medians[species] <- median(bird_lognormals_trunc[species,])
  diffs[species] <- exp(bird_abundances[species,]$meanlog) - means[species]
  diffs_sym[species] <- means[species] - medians[species]
}
sum(means)
mean(diffs)
mean(diffs_sym)




#### Rectified (before summing) Normal
total_dist_recbefore <- c()

for (draw in 1:num_samples){
  total_dist_recbefore[draw] <- sum(bird_lognormals_recbefore[,draw])
}

df_2C_recbefore <- data.frame(sum = total_dist_recbefore)
sd_recbefore <- round(sd(total_dist_recbefore)/1000000000,3) # in billions

p2C_recbefore <- ggplot(df_2C_recbefore, aes(sum)) + geom_density(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 51000000000, color = "black") # + geom_label(aes(x=51000000000, label="51 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recbefore,"billion"), y=5), colour="black")
p2C_recbefore

ggsave("rectifiedbefore_normal.eps")

l2C_recbefore <- ggplot(df_2C_recbefore, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_discrete(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
l2C_recbefore

ggsave("rectifiedbefore_normal_linear.eps")


median(total_dist_recbefore)/1000000000 # in billions
mean(total_dist_recbefore)/1000000000 # in billions




#### Rectified (after summing) Normal
total_dist_recafter <- c()

for (draw in 1:num_samples){
  
  total_dist_recafter[draw] <- sum(bird_lognormals_recafter[,draw])
  
}

median(total_dist_recafter)/1000000000 # in billions
mean(total_dist_recafter)/1000000000 # in billions
hist(total_dist_recafter)

total_dist_recafter <- pmax(0, total_dist_recafter)

median(total_dist_recafter)/1000000000 # in billions
mean(total_dist_recafter)/1000000000 # in billions
hist(total_dist_recafter)




df_2C_recafter <- data.frame(sum = total_dist_recafter)
sd_recafter <- round(sd(total_dist_recafter)/1000000000,3) # in billions

p2C_recafter <- ggplot(df_2C_recafter, aes(sum)) + geom_density(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
p2C_recafter

h2C_recafter <- ggplot(df_2C_recafter, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_log10(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
h2C_recafter

ggsave("rectifiedafter_normal.eps")

l2C_recafter <- ggplot(df_2C_recafter, aes(sum)) + geom_histogram(color = "black", fill="darkgray") + scale_x_discrete(limits = c(2e9, 2e13), breaks = c(1e10, 1e11, 1e12, 1e13), label = c("10,000","100,000","1,000,000","10,000,000")) + xlab("Number of individual birds (millions)") + ylab("") #+ geom_vline(xintercept = 50000000000, color = "black") # + geom_label(aes(x=50000000000, label="50 billion", y=9.5), colour="black") + geom_label(aes(x=5e12, label=paste("sd =",sd_recafter,"billion"), y=5), colour="black")
l2C_recafter

ggsave("rectifiedafter_normal_linear.eps")


median(total_dist_recafter)/1000000000 # in billions
mean(total_dist_recafter)/1000000000 # in billions


p1 <- hist(total_dist_recafter)
plot(p1)
abline(median(total_dist_recafter), add = T)



hist(total_dist_recbefore)
min(total_dist_recbefore)

hist(total_dist_recafter)



# table reporting summary statistics for each alternative solutions
table <- data.frame(method <- NA, median <- NA, mean <- NA, sd <- NA)

table[nrow(table) + 1,] <- c("raw lognormal", round(median(total_dist_raw)/1000000000,3), round(mean(total_dist_raw)/1000000000,3), round(sd(total_dist_raw)/1000000000,3))

table[nrow(table) + 1,] <- c("sorted lognormal", round(median(total_dist_sort)/1000000000,3), round(mean(total_dist_sort)/1000000000,3), round(sd(total_dist_sort)/1000000000,3))

table[nrow(table) + 1,] <- c("shifted lognormal", round(median(total_dist_shift)/1000000000,3), round(mean(total_dist_shift)/1000000000,3), round(sd(total_dist_shift)/1000000000,3))

table[nrow(table) + 1,] <- c("truncated normal", round(median(total_dist_trunc)/1000000000,3), round(mean(total_dist_trunc)/1000000000,3), round(sd(total_dist_trunc)/1000000000,3))

table[nrow(table) + 1,] <- c("rectified normal (before)", round(median(total_dist_recbefore)/1000000000,3), round(mean(total_dist_recbefore)/1000000000,3), round(sd(total_dist_recbefore)/1000000000,3))

table[nrow(table) + 1,] <- c("rectified normal (after)", round(median(total_dist_recafter)/1000000000,3), round(mean(total_dist_recafter)/1000000000,3), round(sd(total_dist_recafter)/1000000000,3))

table <- table[-c(1),]

names(table)[names(table) == 'method....NA'] <- "Method"
names(table)[names(table) == 'median....NA'] <- "median"
names(table)[names(table) == 'mean....NA'] <- "mean"
names(table)[names(table) == 'sd....NA'] <- "sd"

table




bird_lognormals_data <- data.frame(matrix(NA, nrow = num_species, ncol = 8))
colnames(bird_lognormals_data) <- c("CommonName", "ScientificName", "Order", "Family", "MeanLog", "SDLog", "Median", "SD")

for (species in 1:num_species){
  bird_lognormals_data[species,1] <- bird_abundances$CommonName[species]
  bird_lognormals_data[species,2] <- bird_abundances$ScientificName[species]
  bird_lognormals_data[species,3] <- bird_abundances$Order[species]
  bird_lognormals_data[species,4] <- bird_abundances$Family[species]
  
  bird_lognormals_data[species,5] <- bird_abundances$meanlog[species]
  bird_lognormals_data[species,6] <- bird_abundances$sdlog[species]
  bird_lognormals_data[species,7] <- exp(bird_abundances$meanlog[species])
  bird_lognormals_data[species,8] <- exp(bird_abundances$sdlog[species])
}



write.csv(bird_lognormals_data, file = "global_bird_abundance_reconstructed_data.csv")


#### Appendix S1 Fig
sim_mus <- seq(2, 12, 0.1)

sim_df <- data.frame(mus <- c(sim_mus, 10, 12), sigmas <- c(rep(0.2, times = length(sim_mus)), 0.1, 0.4))
colnames(sim_df) <- c("meanlog", "sdlog")

bird_abundance_params <- data.frame(meanlog <- bird_abundances$meanlog, sdlog <- bird_abundances$sdlog)
colnames(bird_abundance_params) <- c("meanlog", "sdlog")

sim_df$id <- "simulation study"
bird_abundance_params$id <- "Callaghan et al. (2021) re-analysis"

df_all <- rbind(sim_df, bird_abundance_params)


ggplot(df_all, aes(x=sdlog, y = meanlog, group = id, colour = id, fill = id)) + geom_point(alpha = 0.5) + xlim(0.05, 2.55)

#ggsave("meanlog_sdlog_scatter.pdf")
