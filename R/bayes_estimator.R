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
  } else if (family = "poisson") {
    fam <- stats::poisson("log")
  }

  vars.response <- all.vars(formula)[1:2]
  vars.predictor <- all.vars(formula)[3:length(all.vars(formula))]

  fit_partialpool <- rstanarm::stan_glmer(formula,
                                          data = data,
                                          family = fam,
                                          ...)

  alphas <- shift_draws(as.matrix(fit_partialpool))
  partialpool <- summary_stats(alphas)
  partialpool <- as.data.frame(partialpool[-nrow(partialpool),])
  colnames(partialpool) <- c("mean", "sd", "Q025", "Q25", "Q50", "Q75", "Q975")

  names <- substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))

  for (i in length(vars.predictor)) {
    partialpool[, vars.predictor[i]] <- unlist(lapply(strsplit(random1, ":"), function(x) x[i+length(vars.predictor)]))
  }

  rownames(partialpool) <- NULL

  partialpool <- suppressWarnings(dplyr::left_join(partialpool, data[c(vars.predictor, "forest")], by = vars.predictor))
  partialpool[, which(!(colnames(partialpool) %in% c(vars.predictor, "forest")))] <- partialpool[, which(!(colnames(partialpool) %in% c(vars.predictor, "forest")))] / partialpool[,"forest"]
  partialpool <- partialpool[, -which(colnames(partialpool) == "forest")]

  return(list(estimate = partialpool,
              model = fit_partialpool))

}
