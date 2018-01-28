#' Function for estimating disturbance rates from TimeSync samples using Bayesian hierarchical modeling
#'
#' @param x A data frame returned from \code{disturbance_summary()}.
#' @param disturbance_col A vector indexing the column storing the number of disturbed plots. Either index value or column name.
#' @param total_col A vector indexing the column storing the number of forested plots. Either index value or column name.
#' @param index_cols A vector indexing the columns storing the hierarchical levels (e.g., years, agents). Either index value or column name.
#' @param prob A vector indicating the quantiles used for summarizing the posterior. Default is 'c(0.025, 0.2, 0.5, 0.8, 0.975)'
#' @param model The model family. Either 'binomial' or 'poisson'.
#' @param trend Wehther a trend should be caluclated.
#' @return A sumary of the joint posterior distribution for each hierarchical level + posterior of parameters (and trend)
#' @export

bayes_estimator <- function (x,
                             disturbance_col,
                             total_col,
                             index_cols,
                             prob = c(0.025, 0.2, 0.5, 0.8, 0.975),
                             model,
                             trend = FALSE) {

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  N <- dim(x)[1]
  K <- as.data.frame(x)[ , total_col]
  y <- as.data.frame(x)[, disturbance_col]
  time <- 1:N

  if (model == "binomial") {
    if (trend) {
      fit <- sampling(stanmodels$bayes_estimator_binomial_trend,
                      data = c("N", "K", "y", "time"),
                      iter = 2000,
                      chains = 4)

      # fit <- stan("src/stan_files/bayes_estimator_binomial_trend.stan",
      #                 data = c("N", "K", "y", "time"),
      #                 iter = 2000,
      #                 chains = 4)

    } else if (trend) {
      fit <- sampling(stanmodels$bayes_estimator_binomial,
                      data = c("N", "K", "y"),
                      iter = 2000,
                      chains = 4)
    }
  } else if (model == "poisson") {
    if (trend) {
      stop("Not implemented yet, sorry...")
    } else if (trend) {
      fit <- sampling(stanmodels$bayes_estimator_poisson,
                      data = c("N", "K", "y"),
                      iter = 2000,
                      chains = 4)
    }
  }

  params <- as.matrix(fit)

  theta <- params[, grep("theta*", colnames(params))]

  if (model == "poisson") {
    theta <- sweep(theta, 2, K, "/")
  }

  estimates <-t(apply(theta, 2, function(z) {c(mean(z), sd(z), quantile(z, prob))}))
  estimates <- as.data.frame(estimates)
  colnames(estimates) <- c("mean", "sd", paste0("Q", prob))
  for (i in index_cols) estimates[, i] <- x[, i]
  rownames(estimates) <- NULL

  posterior <- data.table::melt(theta)
  for (i in index_cols) posterior[, i] <- rep(x[, i][[1]], each = length(unique(posterior$iterations)))
  posterior$parameters <- NULL

  if (trend) {

    trend_pred <- params[, grep("trend_pred*", colnames(params))]
    trend <- params[, grep("trend$", colnames(params))]

    trend_estimates <-t(apply(trend_pred, 2, function(z) {c(mean(z), sd(z), quantile(z, prob))}))
    trend_estimates <- as.data.frame(trend_estimates)
    colnames(trend_estimates) <- c("mean", "sd", paste0("Q", prob))
    for (i in index_cols) trend_estimates[, i] <- x[, i]
    rownames(trend_estimates) <- NULL

    trend_posterior <- data.table::melt(trend_pred)
    for (i in index_cols) trend_posterior[, i] <- rep(x[, i][[1]], each = length(unique(trend_posterior$iterations)))
    trend_posterior$parameters <- NULL

    return(list(estimate = estimates,
                posterior = posterior,
                model = fit,
                trend_estimate = trend_estimates,
                trend_posterior = trend_posterior,
                trend = trend))

  } else if (trend) {

    return(list(estimate = estimates,
                posterior = posterior,
                model = fit))

  }

}

