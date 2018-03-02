
// Stan model for partial pooling using a hierarchical beta model (partially-pooled binary trials) with linear trend component

data {
  int<lower=0> N;                       // items
  int<lower=0> K[N];                    // initial forested plots
  int<lower=0> y[N];                    // initial disturbance plots
  int<lower=0> time[N];                 // time signal for trend prediction
  int<lower=0> length_pred;             // length of prediction
  int<lower=0> time_pred[length_pred];  // time signal for predictions
}

transformed data {

  // Dervide statistics for Bayesian p-value calculations

  real min_y;   // minimum successes
  real max_y;   // maximum successes
  real mean_y;  // sample mean successes
  real sd_y;    // sample std dev successes

  min_y = min(y);
  max_y = max(y);
  mean_y = mean(to_vector(y));
  sd_y = sd(to_vector(y));
}

parameters {
  real<lower=0> sigma;          // population sd of disturbance log-odds
  real trend;                   // trend in disturbance rate
  real intercept;               // intercept
  vector[N] alpha_std;          // disturbance log-odds
}

transformed parameters {
  vector[N] mu;
  for (i in 1:N)
    mu[i] = intercept + trend * time[i];
}

model {
  trend ~ student_t(5, 0, 0.1);                  // hyperprior
  intercept ~ normal(-6, 1);                      // hyperprior
  sigma ~ normal(0, 1);                          // hyperprior
  alpha_std ~ normal(0, 1);                      // prior

  y ~ binomial_logit(K, mu + sigma * alpha_std); // likelihood
}

generated quantities {

  vector[N] theta;  // chance of success
  vector[length_pred] trend_pred;  // trend predictions
  vector[N] log_lik; // pointwise log-likelihood

  int<lower=0> y_rep[N];      // replications for existing items
  int<lower=0> y_pop_rep[N];  // replications for simulated items

  real<lower=0> min_y_rep;   // posterior predictive min replicated successes
  real<lower=0> max_y_rep;   // posterior predictive max replicated successes
  real<lower=0> mean_y_rep;  // posterior predictive sample mean replicated successes
  real<lower=0> sd_y_rep;    // posterior predictive sample std dev replicated successes

  int<lower=0, upper=1> p_min;  // posterior predictive p-values
  int<lower=0, upper=1> p_max;
  int<lower=0, upper=1> p_mean;
  int<lower=0, upper=1> p_sd;

  theta = inv_logit(mu + sigma * alpha_std);

  for (i in 1:length_pred)
    trend_pred[i] = inv_logit(intercept + trend * time_pred[i]);

  for (n in 1:N)
    log_lik[n] = binomial_logit_lpmf(y[n] | K[n], mu + sigma * alpha_std);

  for (n in 1:N)
    y_rep[n] = binomial_rng(K[n], theta[n]);

  for (n in 1:N)
    y_pop_rep[n] = binomial_rng(K[n], inv_logit(normal_rng(mu[n], sigma)));

  min_y_rep = min(y_rep);
  max_y_rep = max(y_rep);
  mean_y_rep = mean(to_vector(y_rep));
  sd_y_rep = sd(to_vector(y_rep));

  p_min = (min_y_rep >= min_y);
  p_max = (max_y_rep >= max_y);
  p_mean = (mean_y_rep >= mean_y);
  p_sd = (sd_y_rep >= sd_y);

}
