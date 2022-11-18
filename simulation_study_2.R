
#### Talis et al
# Difficulties in summing distributions for abundance and potential solutions

# simulation study 2

library(compositions)
library(tidyverse)
library(dfrtopics)
# devtools::install_github("agoldst/dfrtopics")
library(abind)
library(ggplot2)


# n = populations and i = "posterior" samples from each population
# mu is mean of logged abundance and sd is standard deviation of logged abundance
mu <- 100
sd <- .2
n <- 100
i <- 1000

myfunc <- function(phi = 0) {
  x <- matrix(NA, nrow = n, ncol = i)
  v <- matrix(sd * sd * phi, nrow = n, ncol = n)
  diag(v) <- sd^2
  x <- matrix(compositions::rlnorm.rplus(i, meanlog = rep(log(mu), n), varlog = v), 
              nrow = n, ncol = i, byrow = T)
  return(median(colSums(x)) / sum(apply(x, 1, median)))
}

phi_vec <- seq(0, 1, .05)
bias_matrix <- matrix(NA, nrow = 500, ncol = length(phi_vec))

for (k in 1:(length(phi_vec))) {
  if (k != length(phi_vec)) {
    for (m in 1:500) {
      bias_matrix[m, k] <- myfunc(phi = phi_vec[k])
    }
  } else {
    for (m in 1:500) {
      pop <- rlnorm(i, meanlog = log(mu), sdlog = sd)
      x <- matrix(pop, nrow = n, ncol = i, byrow = TRUE)
      bias_matrix[m, k] <- median(colSums(x)) / sum(apply(x, 1, median))
    }
  }    
}

bias_df <- dfrtopics::gather_matrix(bias_matrix) %>%
  mutate(phi = rep(phi_vec, 500))
ggplot(bias_df, aes(x = as.factor(phi), y = value)) +
  geom_boxplot(fill = '#A4A4A4', color = "black") +
  labs(x = "correlation coefficient", y = "median(sums) / sum(medians)") +
  theme_minimal()


save(bias_df, file = "simulation_study_2.RData")
