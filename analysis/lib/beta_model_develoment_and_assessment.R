
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("R/disturbance_summary.R")
source("R/bayes_estimator_year.R")

dat <- read.csv("analysis/data/austria_vertex_07_07_17.csv")

dat_summary <- disturbance_summary(dat, by.year = TRUE, by.agent = FALSE)

fit <- bayes_estimator_year(dat_summary, c(0.025, 0.2, 0.25, 0.5, 0.75, 0.8, 0.975))

ggplot(fit$estimate, aes(x = year, y = Q0.5 * 100)) +
  geom_ribbon(aes(ymin = Q0.2 * 100, ymax = Q0.8 * 100), alpha = 0.3) +
  geom_line()
