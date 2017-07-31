#' Funktion for comparing models
#'
#' @param newdata Newdata.
#' @param models Models.
#' @export

compare_models <- function(models, newdata = NULL, density = FALSE) {

  log_p_new_mats <- lapply(models, rstanarm::log_lik, newdata = newdata)

  if (density) {

    log_sum_exp <- function(u) {
      max_u <- max(u)
      max_u + log(sum(exp(u - max_u)))
    }

    log_mean_exp <- function(u) {
      M <- length(u)
      -log(M) + log_sum_exp(u)
    }

    new_lps <- lapply(log_p_new_mats, function(x) apply(x, 2, log_mean_exp))
    new_lps_sums <- sapply(new_lps, sum)

    return(new_lps_sums)

  } else {

    log_p_new <- sapply(log_p_new_mats, rowSums)
    mean_log_p_new <- colMeans(log_p_new)

    return(mean_log_p_new)
  }

}
