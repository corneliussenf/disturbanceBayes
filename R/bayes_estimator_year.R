
bayes_estimator_year <- function (x, p, model, omission_rate = NULL) {

  require(rstan)

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  N <- dim(x)[1]
  K <- x$forest
  y <- x$disturbance

  if (model == "binomial") {
    fit <- rstan::stan("Stan/bayes_estimator_binomial.stan",
                       data = c("N", "K", "y"),
                       iter = 2000,
                       chains = 4)
  } else if (model == "poisson") {
    fit <- rstan::stan("Stan/bayes_estimator_poisson.stan",
                       data = c("N", "K", "y"),
                       iter = 2000,
                       chains = 4)
  }

  theta <- as.matrix(fit)
  theta <- theta[, grep("theta*", colnames(theta))]

  if (model == "poisson") {
    theta <- sweep(theta, 2, K, "/")
  }

  estimates <-t(apply(theta, 2, function(z) {c(mean(z), sd(z), quantile(z, p))}))
  estimates <- as.data.frame(estimates)
  colnames(estimates) <- c("mean", "sd", paste0("Q", p))
  estimates$year <- x$year
  rownames(estimates) <- NULL

  posterior <- data.table::melt(theta)
  posterior$year <- rep(x$year, each = length(unique(posterior$iterations)))
  posterior$parameters <- NULL

  return(list(estimate = estimates,
              posterior = posterior,
              model = fit))

}

