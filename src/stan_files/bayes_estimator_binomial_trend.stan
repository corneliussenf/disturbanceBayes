
// Stan model for partial pooling using a hierarchical beta model (partially-pooled binary trials) with linear trend component

data {
  int<lower=0> N;           // items
  int<lower=0> K[N];        // initial forested plots
  int<lower=0> y[N];        // initial disturbance plots
  int<lower=0> time[N];
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
  trend ~ normal(0, 0.5);                          // hyperprior
  intercept ~ normal(1, 1);                      // hyperprior
  sigma ~ normal(0, 1);                          // hyperprior
  alpha_std ~ normal(0, 1);                      // prior

  y ~ binomial_logit(K, mu + sigma * alpha_std); // likelihood
}

generated quantities {
  vector[N] theta;  // chance of success
  vector[N] trend_pred;  // trend predictions
  vector[N] log_lik; // pointwise log-likelihood
  theta = inv_logit(mu + sigma * alpha_std);
  for (i in 1:N)
    trend_pred[i] = inv_logit(intercept + trend * time[i]);
  for (i in 1:N)
    log_lik[i] = binomial_logit_lpmf(y[i] | K[i], mu + sigma * alpha_std);
}
