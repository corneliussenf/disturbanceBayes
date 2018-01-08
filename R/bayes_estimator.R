#' Function for estimating disturbance rates from TimeSync samples using Bayesian hierarchical modeling
#'
#' @param formula A model formula.
#' @param data A data.frame containing the data.
#' @param family The model family. Currently implemented are 'binomial' and 'poisson'.
#' @param n_cores The number of cores to be used for sampling the posterior.
#' @return A sumary of the joint posterior distribution.
#' @export

bayes_estimator <- function(formula,
                            data,
                            family,
                            n_cores = 4,
                            weights = NULL,
                            ...) {

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = n_cores)

  invlogit <- function(x) 1/(1 + exp(-x))

  shift_draws <- function(draws) {
    sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
  }

  summary_stats <- function(posterior, f = family) {
    if (f == "binomial") {
      x <- invlogit(posterior)
    } else if (f == "poisson") {
      x <- exp(posterior)
    }

    t(apply(x, 2,  function(x) {c(mean = mean(x),
                     sd = sd(x),
                     Q025 = quantile(x, 0.025),
                     Q25 = quantile(x, 0.25),
                     Q50 = quantile(x, 0.5),
                     Q75 = quantile(x, 0.75),
                     Q975 = quantile(x, 0.975))}))
  }

  if (family == "binomial") {
    fam <- stats::binomial("logit")
  } else if (family == "poisson") {
    fam <- stats::poisson("log")
  }

  if (family == "binomial") {
    vars.response <- all.vars(formula)[1:2]
    vars.predictor <- all.vars(formula)[3:length(all.vars(formula))]
  } else if (family == "poisson") {
    vars.response <- all.vars(formula)[1]
    vars.predictor <- all.vars(formula)[2:length(all.vars(formula))]
  }

  fit <- rstanarm::stan_glmer(formula,
                              data = data,
                              family = fam,
                              weights = weights,
                              ...)

  alphas <- shift_draws(as.matrix(fit))
  estimates <- summary_stats(alphas)
  estimates <- as.data.frame(estimates[-nrow(estimates),])
  colnames(estimates) <- c("mean", "sd", "Q025", "Q25", "Q50", "Q75", "Q975")

  posterior <- data.table::melt(alphas[,-ncol(alphas)])
  if (family == "binomial") {
    posterior$value <- invlogit(posterior$value)
  } else if (family == "poisson") {
    posterior$value <- exp(posterior$value)
  }

  names1 <- substring(rownames(estimates), 15, (nchar(rownames(estimates)) - 1))
  names2 <- substring(as.character(posterior$parameters), 15, (nchar(as.character(posterior$parameters)) - 1))

  for (i in 1:length(vars.predictor)) {
    labels1 <- unlist(lapply(strsplit(names1, ":"), function(x) x[i + length(vars.predictor)]))
    estimates[, vars.predictor[i]] <- labels1

    labels2 <- unlist(lapply(strsplit(names2, ":"), function(x) x[i + length(vars.predictor)]))
    posterior[, vars.predictor[i]] <- labels2
  }

  rownames(estimates) <- NULL
  posterior <- posterior[, -which(colnames(posterior) == "parameters")]

  if (family == "poisson") {
    estimates <- suppressWarnings(dplyr::left_join(estimates, data[c(vars.predictor, "forest")], by = vars.predictor))
    estimates[, which(!(colnames(estimates) %in% c(vars.predictor, "forest")))] <- estimates[, which(!(colnames(estimates) %in% c(vars.predictor, "forest")))] / estimates[,"forest"]
    estimates <- estimates[, -which(colnames(estimates) == "forest")]

    posterior <- suppressWarnings(dplyr::left_join(posterior, data[c(vars.predictor, "forest")], by = vars.predictor))
    posterior[, "value"] <- posterior[, "value"] / posterior[,"forest"]
    posterior <- posterior[, -which(colnames(posterior) == "forest")]
  }

  return(list(estimate = estimates,
              posterior = posterior,
              model = fit))

}
