#' Function for summarizing disturbances from a TimeSync sample
#'
#' @param dat Data exported from TimeSync.
#' @param by.agent Logical. Whether disturbances should be aggregated by agent.
#' @param grouping.vars Data.frame. Additional grouping variables stored in a data frame, where each row corresponds to one plot. The data frame must have a column 'plotid'.
#' @param change.agents Change agents to be considered. Standard is 'Harvest', 'Wind', and 'Decline'.
#' @return A sumary of the disturbances per year (and agents/grouping variables).
#' @export
#'
disturbance_summary <- function(dat, by.agent = FALSE, grouping.vars = NULL, change.agents = c("Harvest", "Wind", "Decline")) {

  dat_processed <- dplyr::mutate(dat, agent = lead(change_process))
  dat_processed <- dplyr::filter(dat_processed, agent %in% change.agents & dominant_landuse == "Forest")
  dat_processed <- dplyr::select(dat_processed, -change_process)

  if (is.null(grouping.vars)) {
    if (by.agent) {
      dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, "image_year", "agent"), disturbance = length(agent))
    } else {
      dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, "image_year"), disturbance = length(agent))
    }
  } else {
    dat_processed <- dplyr::left_join(dat_processed, grouping.vars, by = "plotid")
    grouping.vars.names <- names(grouping.vars)[-which(names(grouping.vars) == "plotid")]
    if (by.agent) {
      dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, .dots = c("image_year", lapply(grouping.vars.names, function(x) x ), "agent")), disturbance = length(agent))
    } else {
      dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, .dots = c("image_year", lapply(grouping.vars.names, function(x) x ))), disturbance = length(agent))
    }
  }


  forest <- sum(dplyr::summarize(dplyr::group_by(dat, plotid), forest = sum(dominant_landuse == "Forest") > 0)$forest)
  dat_processed <- dplyr::mutate(dat_processed, forest = forest, year = image_year + 1)
  dat_processed <- dplyr::ungroup(dat_processed)
  dat_processed <- dplyr::select(dat_processed, -image_year)

  return(dat_processed)
}
