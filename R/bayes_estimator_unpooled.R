#' Function for estimating disturbance rates from TimeSync samples using Bayesian modeling
#'
#' @param formula A model formula.
#' @param data A data.frame containing the data.
#' @param family The model family. Currently implemented are 'binomial' and 'poisson'.
#' @param n_cores The number of cores to be used for sampling the posterior.
#' @return A sumary of the joint posterior distribution.
#' @export

bayes_estimator_unpooled <- function(formula,
                                     data,
                                     family,
                                     n_cores = 4,
                                     ...) {

  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = n_cores)

  invlogit <- function(x) 1/(1 + exp(-x))

  shift_draws <- function(draws) {
    out <- sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
    cbind(draws[, 1], out)
  }

  summary_stats <- function(posterior, f = family) {
    if (f == "binomial") {
      x <- invlogit(posterior)
    } else if (f == "poisson") {
      x <- exp(posterior)
    }

    t(apply(x, 2, function(x) {c(mean = mean(x),
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

  fit <- rstanarm::stan_glm(formula,
                            data = data,
                            family = fam,
                            ...)

  alphas <- shift_draws(as.matrix(fit))
  estimate <- summary_stats(alphas)
  estimate <-as.data.frame(estimate)
  colnames(estimate) <- c("mean", "sd", "Q025", "Q25", "Q50", "Q75", "Q975")

  posterior <- data.table::melt(alphas)
  if (family == "binomial") {
    posterior$value <- invlogit(posterior$value)
  } else if (family == "poisson") {
    posterior$value <- exp(posterior$value)
  }
  posterior$iterations <- posterior$Var1

  if (length(vars.predictor) == 1) {
    levels <- gsub(vars.predictor, "", rownames(estimate))
    names <- unique(data[,vars.predictor])
    levels <- c(names[which(!(names %in% levels))], levels[-1])

    estimate[, vars.predictor] <- levels
    rownames(estimate) <- NULL

    posterior[, vars.predictor] <- levels
    posterior <- posterior[, -c(1:2)]

  } else {

    for (i in 1:length(vars.predictor)) {
      names <- unlist(lapply(strsplit(rownames(estimate), ":"), function(x) x[i]))
      estimate[, vars.predictor[i]] <- gsub(vars.predictor[i], "", names)
      posterior[, vars.predictor[i]] <- gsub(vars.predictor[i], "", names)
    }

    estimate <- estimate[-1, ]
    rownames(estimate) <- NULL

    posterior <- posterior[, -c(1:2)]

  }

  if (family == "poisson") {
    estimate <- suppressWarnings(dplyr::left_join(estimate, data[c(vars.predictor, "forest")], by = vars.predictor))
    estimate[, which(!(colnames(estimate) %in% c(vars.predictor, "forest")))] <- estimate[, which(!(colnames(estimate) %in% c(vars.predictor, "forest")))] / estimate[,"forest"]
    estimate <- estimate[, -which(colnames(estimate) == "forest")]

    posterior <- suppressWarnings(dplyr::left_join(posterior, data[c(vars.predictor, "forest")], by = vars.predictor))
    posterior[, "value"] <- posterior[, "value"] / posterior[,"forest"]
    posterior <- posterior[, -which(colnames(posterior) == "forest")]
  }

  return(list(estimate = estimate,
              posterior = posterior,
              model = fit))

}
