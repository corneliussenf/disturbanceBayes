
// Stan model for partial pooling using a hierarchical beta model (partially-pooled binary trials)

data {
  int<lower=0> N;           // items
  int<lower=0> K[N];        // initial forested plots
  int<lower=0> y_obs[N];    // initial disturbance plots
}

parameters {
  real mu;                       // population mean of disturbance log
  real<lower=0> sigma;           // population sd of disturbance log
  vector[N] beta_std;            // disturbance log-odds
}

model {
  mu ~ normal(1, 1);                             // hyperprior
  sigma ~ normal(0, 1);                          // hyperprior
  beta_std ~ normal(0, 1);                       // prior

  y ~ poisson_log(mu + sigma * beta_std);  // likelihood
}

generated quantities {

  vector[N] theta;  // chance of success
  vector[N] log_lik; // pointwise log-likelihood

  theta = exp(mu + sigma * beta_std);
  for (i in 1:N)
    log_lik[i] = poisson_log_lpmf(y[i] | mu + sigma * beta_std);

}
