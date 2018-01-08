
library(tidyverse)
#devtools::install_github("corneliussenf/disturbanceBayes")
library(disturbanceBayes)

dat <- read.csv("../../TimeSync/data/timeSync/switzerland/ts_vertex_switzerland_141217.csv")
#dat <- read.csv("../../TimeSync/data/timeSync/slovakia/ts_vertex_17102017.csv")

dat_summary <- disturbance_summary(dat, by.year = TRUE, by.agent = FALSE)

fit_poisson <- bayes_estimator_year(dat_summary,
                                    c(0.025, 0.2, 0.25, 0.5, 0.75, 0.8, 0.975),
                                    "poisson")

fit_binomial <- bayes_estimator_year(dat_summary,
                                     c(0.025, 0.2, 0.25, 0.5, 0.75, 0.8, 0.975),
                                     "binomial")

ggplot() +
  geom_ribbon(data = fit_poisson$estimate,
              aes(x = as.integer(year),
                  ymin = Q0.25 * 100,
                  ymax = Q0.75 * 100), alpha = 0.4, fill = "blue") +
  geom_ribbon(data = fit_binomial$estimate,
              aes(x = as.integer(year),
                  ymin = Q0.25 * 100,
                  ymax = Q0.75 * 100), alpha = 0.4, fill = "red") +
  geom_line(data = fit_poisson$estimate, aes(x = as.integer(year), y = Q0.5 * 100), col = "blue") +
  geom_line(data = fit_binomial$estimate, aes(x = as.integer(year), y = Q0.5 * 100), col = "red") +
  labs(x = "Year", y = bquote("Disturbance rate [% "*yr^-1*"]")) +
  ggthemes::theme_few()

log_lik_poisson <- loo::extract_log_lik(fit_poisson$model, parameter_name = "log_lik")
log_lik_poisson <- apply(log_lik_poisson, 1, sum)

log_lik_binomial <- loo::extract_log_lik(fit_binomial$model, parameter_name = "log_lik")
log_lik_binomial <- apply(log_lik_binomial, 1, sum)

#a <- exp(as.brob(log_lik_binomial)) / (exp(as.brob(log_lik_poisson)) + exp(as.brob(log_lik_binomial)))
a <- log_lik_binomial / (log_lik_poisson + log_lik_binomial)

hist(a)

which <- vector("integer", length(a))
for (i in 1:length(a)) which[i] <- sample(c(1, 2), size = 1, replace = TRUE, prob = c(as.double(a)[i], 1 - as.double(a)[i]))

dd <- data.frame(binomial = fit_binomial$posterior$value,
                 poisson = fit_poisson$posterior$value,
                 which = which)
dd <- apply(dd, 1, function(x) x[x[3]])

ddd <- data.frame(iteration = fit_binomial$posterior$iterations,
                  year = fit_binomial$posterior$year,
                  value = dd)

ddd_summary <- ddd %>%
  group_by(year) %>%
  summarize(estimate = median(value),
            lower = quantile(value, 0.2),
            upper = quantile(value, 0.8))

ggplot() +
  geom_line(data = fit_poisson$estimate, aes(x = as.integer(year), y = Q0.5 * 100), col = "blue") +
  geom_line(data = fit_binomial$estimate, aes(x = as.integer(year), y = Q0.5 * 100), col = "red") +
  geom_line(data = ddd_summary, aes(x = year, y = estimate * 100))



