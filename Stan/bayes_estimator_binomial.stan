
// Stan model for partial pooling using a hierarchical beta model (partially-pooled binary trials)

data {
  int<lower=0> N;           // items
  int<lower=0> K[N];        // initial forested plots
  int<lower=0> y[N];        // initial disturbance plots
  // real omission_rate;
}

// transformed data {
//   int omission[N];
//
//   for (i in 1:N)
//     omission[i] = binomial_rng(y[i], omission_rate);
//
// }

parameters {
  real mu;                       // population mean of disturbance log-odds
  real<lower=0> sigma;           // population sd of disturbance log-odds
  vector[N] alpha_std;           // disturbance log-odds
}

model {

  mu ~ normal(1, 1);                             // hyperprior
  sigma ~ normal(0, 1);                           // hyperprior
  alpha_std ~ normal(0, 1);                      // prior

  y ~ binomial_logit(K, mu + sigma * alpha_std);  // likelihood
}

generated quantities {

  vector[N] theta;  // chance of success
  vector[N] log_lik; // pointwise log-likelihood

  theta = inv_logit(mu + sigma * alpha_std);
  for (i in 1:N)
    log_lik[i] = binomial_logit_lpmf(y[i] | K[i], mu + sigma * alpha_std);

}
