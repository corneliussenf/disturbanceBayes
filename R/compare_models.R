#' Add together two numbers
#'
#' @param formula A model formula.
#' @param data A data.frame containing the data.
#' @param family The model family. Currently implemented are 'binomial' and 'poisson'.
#' @param summary_probs Vector containing the quantiles for summarizing the joint posterior distribution.
#' @return A sumary of the joint posterior distribution.
#' @export

compare_models <- function(newdata, models) {
  log_p_new_mats <- lapply(models, log_lik, newdata = newdata)
  log_p_new <- sapply(log_p_new_mats, rowSums)
  M <- nrow(log_p_new)
  mean_log_p_new <- colMeans(log_p_new)
  round(sort(mean_log_p_new, decreasing = TRUE), digits = 1)
}

