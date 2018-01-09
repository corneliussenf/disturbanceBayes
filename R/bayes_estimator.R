#' Function for estimating disturbance rates from TimeSync samples using Bayesian hierarchical modeling
#'
#' @param x A data frame returned from \code{disturbance_summary()}.
#' @param p A vector indicating the quantiles used for summarizing the posterior.
#' @param index_cols A vector indexing the columns storing the hierarchical levels (e.g., years). Either index value or column name.
#' @param model The model family. Either 'binomial' or 'poisson'.
#' @return A sumary of the joint posterior distribution for each hierarchical level.
#' @export

bayes_estimator <- function (x, p, index_cols, model) {

  require(rstan)

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  N <- dim(x)[1]
  K <- x$forest
  y <- x$disturbance

  if (model == "binomial") {
    fit <- rstan::sampling(stanmodels$bayes_estimator_binomial,
                       data = c("N", "K", "y"),
                       iter = 2000,
                       chains = 4)
  } else if (model == "poisson") {
    fit <- rstan::sampling(stanmodels$bayes_estimator_poisson,
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
  for (i in index_cols) estimates[, i] <- x[, i]
  rownames(estimates) <- NULL

  posterior <- data.table::melt(theta)
  posterior$year <- rep(x$year, each = length(unique(posterior$iterations)))
  posterior$parameters <- NULL

  return(list(estimate = estimates,
              posterior = posterior,
              model = fit))

}

