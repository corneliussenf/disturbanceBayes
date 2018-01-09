
library(tidyverse)

source("R/disturbance_summary.R")
source("R/bayes_estimator2.R")

dat <- read.csv("../../TimeSync/data/timesync.csv")

dat_summary <- dat %>%
  split(.$project_code) %>%
  map(~ disturbance_summary(., by.agent = FALSE)) %>%
  bind_rows(.id = "country") %>%
  mutate(country = factor(country, labels = c("Austria", "Germany", "Poland", "Slovakia", "Switzerland")))

fit <- bayes_estimator_year(dat_summary, p = c(0.025, 0.2, 0.25, 0.5, 0.75, 0.8, 0.975), model = "binomial")

fit$estimate$year <- dat_summary$year
fit$estimate$country <- dat_summary$country

ggplot(fit$estimate, aes(x = year, y = Q0.5 * 100)) +
  geom_ribbon(aes(ymin = Q0.2 * 100, ymax = Q0.8 * 100), alpha = 0.3) +
  geom_line() +
  facet_wrap(~country)
