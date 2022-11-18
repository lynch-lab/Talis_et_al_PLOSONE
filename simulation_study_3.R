
#### Talis et al
# Difficulties in summing distributions for abundance and potential solutions

# simulation study 3

library(compositions)
library(tidyverse)
library(dfrtopics)
# devtools::install_github("agoldst/dfrtopics")
library(abind)
library(patchwork)

mu1 <- 10
mu2 <- 12
sd1 <- .1
sd2 <- .4
pop <- 25
sample <- 5000
sim <- 100
x1 <- array(NA, dim = c(pop, sample, sim))
x2 <- array(NA, dim = c(pop, sample, sim))
for (i in 1:pop) {
  for (j in 1:sim) { 
    x1[i,,j] <- rlnorm(sample, meanlog = mu1, sdlog = sd1)
    x2[i,,j] <- rlnorm(sample, meanlog = mu2, sdlog = sd2)
  }
}
bias_matrix <- array(NA, dim = c(pop + 1, pop + 1))
for (i in 0:pop) {
  for (j in 0:pop) {
    xc <- abind(x1[0:i,,, drop = F], x2[0:j,,, drop = F], along = 1)
    x_sum_medians <- apply(xc, c(1, 3), median)
    x_sum_medians <- apply(x_sum_medians, 2, sum)
    x_median_sums <- apply(xc, c(2, 3), sum)
    x_median_sums <- apply(x_median_sums, 2, median)
    bias_matrix[i + 1, j + 1] <- mean(x_median_sums / x_sum_medians)
  }
}  


bias_df <- dfrtopics::gather_matrix(bias_matrix) %>%
  mutate(population1 = row_key - 1, population2 = col_key - 1) %>%
  mutate(bias = round(100 * (value - 1), 2)) %>%
  dplyr::select(population1, population2, bias) %>%
  mutate(bias_label = sprintf("%.1f", bias))
bias_df$bias_label[1] <- ""
mu_LN <- c(mu1, mu2)
sd_LN <- c(sd1, sd2)
bias_x_label <- bquote("number of populations with" ~ mu[1] * 
                         "=" * .(mu_LN[1]) * "," ~ sigma[1] * "=" * .(sd_LN[1]))
bias_y_label <- bquote("number of populations with" ~ mu[2] * 
                         "=" * .(mu_LN[2]) * "," ~ sigma[2] * "=" * .(sd_LN[2]))
g2 <- ggplot(data = bias_df, aes(x = population1, y = population2, fill = bias, label = bias_label)) +
  geom_tile() + 
  labs(x = bias_x_label, y = bias_y_label, fill = "Percent Bias") +
  scale_fill_gradient2(low = scales::muted("red"), high = scales::muted("blue"), 
                       limits = c(0, max(bias_df$bias))) +
  geom_text(size = 2) +
  theme_classic() +
  scale_x_discrete(limits = seq(0, pop, 1)) +
  scale_y_discrete(limits = seq(0, pop, 1)) +
  theme(axis.title.x = element_text(margin = margin(t = 10))) +
  theme(axis.title.y = element_text(margin = margin(r = 10)))
g2

save(bias_df, file = "simulation_study_3.RData")
