#' Funktion for comparing models
#'
#' @param newdata Newdata.
#' @param models Models.
#' @export

compare_models <- function(newdata, models) {
  log_p_new_mats <- lapply(models, rstanarm::log_lik, newdata = newdata)
  log_p_new <- sapply(log_p_new_mats, rowSums)
  M <- nrow(log_p_new)
  mean_log_p_new <- colMeans(log_p_new)
  round(sort(mean_log_p_new, decreasing = TRUE), digits = 1)
}

