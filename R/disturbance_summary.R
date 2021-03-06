#' Function for summarizing disturbances from a TimeSync sample
#'
#' @param dat Data exported from TimeSync.
#' @param by.year TBD...
#' @param by.agent Logical. Whether disturbances should be aggregated by agent.
#' @param standreplacing Logical. Whether disturbances should be aggregated by stand-replacing/non-stand-replacing
#' @param grouping.vars Data.frame. Additional grouping variables stored in a data frame, where each row corresponds to one plot. The data frame must have a column 'plotid'.
#' @param sub.agents TBD...
#' @param change.agents Change agents to be considered. Standard is 'Harvest', 'Wind', and 'Decline'.
#' @param agent.regroup TBD...
#' @param interpreter TBD...
#' @param remove.1985.decline TBD...
#' @return A sumary of the disturbances per year (and agents/grouping variables).
#' @export

disturbance_summary <- function(dat,
                                by.year = TRUE,
                                by.agent = FALSE,
                                standreplacing = FALSE,
                                grouping.vars = NULL,
                                sub.agents = FALSE,
                                change.agents = c("Harvest", "Wind", "Decline", "Hydrology", "Debris", "Other"),
                                agent.regroup = NULL,
                                interpreter = NULL,
                                remove.1985.decline = TRUE) {

  dat <- dplyr::as_tibble(dat)

  if (!is.null(interpreter)) dat <- dat[dat$interpreter == interpreter, ]

  if (sub.agents) {
    dat$change_process <- as.character(dat$change_process)
    agent_new <- unlist(lapply(strsplit(as.character(dat$change_process_notes), "|", fixed = TRUE), function(x) x[1]))
    dat[!is.na(agent_new) & dat$change_process == "Harvest", "change_process"] <- agent_new[!is.na(agent_new) & dat$change_process == "Harvest"]

    still_harvest <- dat[dat$change_process == "Harvest", "plotid"]
    still_harvest <- na.omit(still_harvest)
    if (nrow(still_harvest) > 0) print(paste0("Warning: Detected 'harvests' without sub-agent, which were automatically set to 'clearcut': ",
                       paste(still_harvest[[1]], collapse = ", ")))
    dat[which(dat$change_process == "Harvest"), "change_process"] <- "clearcut"

    change.agents <- c(change.agents, unique(agent_new))
    change.agents <- change.agents[!is.na(change.agents)]
  }

  dat_processed <- dplyr::mutate(dat,
                                 agent = dplyr::lead(change_process),
                                 post_disturbance_lc = dplyr::lead(dominant_landcover))
  dat_processed <- dplyr::filter(dat_processed, agent %in% change.agents & dominant_landuse == "Forest")
  dat_processed$agent <- as.factor(as.character(dat_processed$agent))
  if (remove.1985.decline == TRUE) dat_processed <- dplyr::filter(dat_processed, !(agent == "Decline" & image_year == 1984))

  if (standreplacing) {
    dat_processed$agent <- as.character(dat_processed$agent)
    dat_processed[which(!is.na(dat_processed$agent) & dat_processed$post_disturbance_lc == "Trees"), "agent"] <- "Non-standreplacing"
    dat_processed[which(dat_processed$agent != "Stable" & dat_processed$post_disturbance_lc == "Non-tree"), "agent"] <- "Standreplacing"
    dat_processed$agent <- as.factor(dat_processed$agent)
  }

  if (!is.null(agent.regroup)) {
    dat_processed <- dplyr::left_join(dat_processed, agent.regroup, by = "agent")
    dat_processed <- dplyr::mutate(dat_processed, "agent" = new)
    dat_processed <- dplyr::select(dat_processed, -new)
  }

  if (is.null(grouping.vars)) {
    dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, "image_year", "agent"), disturbance = length(agent))
    forest <- sum(dplyr::summarize(dplyr::group_by(dat, plotid),
                                   forest = sum(dominant_landuse == "Forest", na.rm = TRUE) > 0)$forest)
    dat_processed <- dplyr::mutate(dat_processed, forest = forest)
  } else {

    dat_processed <- dplyr::left_join(dat_processed, grouping.vars, by = "plotid")
    grouping.vars.names <- names(grouping.vars)[-which(names(grouping.vars) == "plotid")]

    dat <- dplyr::left_join(dat, grouping.vars, by = "plotid")
    forest <- dplyr::summarize(dplyr::group_by_(dat, .dots = unlist(c("plotid", lapply(grouping.vars.names, function(x) x )))),
                                   forest = sum(dominant_landuse == "Forest") > 0)
    forest <- dplyr::summarize(dplyr::group_by_(forest, .dots = unlist(c(lapply(grouping.vars.names, function(x) x )))),
                               forest = sum(forest, na.rm = TRUE))
    forest <- na.omit(forest)
    dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, .dots = c("image_year", lapply(grouping.vars.names, function(x) x ), "agent")), disturbance = length(agent))
    dat_processed <- dplyr::left_join(dat_processed, forest, by = grouping.vars.names)
  }

  dat_processed <- dplyr::mutate(dat_processed, year = image_year + 1)
  dat_processed <- dplyr::ungroup(dat_processed)
  dat_processed <- dplyr::select(dat_processed, -image_year)
  dat_processed$agent <- tolower(dat_processed$agent)

  if (by.agent == FALSE & standreplacing == FALSE) {
    g <- names(dat_processed)[-which(names(dat_processed) %in% c("forest", "disturbance", "agent"))]
    dat_processed <- dplyr::summarize(dplyr::group_by_(dat_processed, .dots = lapply(g, function(x) x)),
                                      disturbance = sum(disturbance),
                                      forest = unique(forest))
  }

  if (by.year == FALSE) {
    g <- names(dat_processed)[-which(names(dat_processed) %in% c("forest", "disturbance", "year"))]
    dat_processed <- dplyr::summarize(dplyr::group_by_(dat_processed, .dots = lapply(g, function(x) x)),
                                      disturbance = sum(disturbance),
                                      forest = unique(forest))
  }

  dat_processed <- dplyr::ungroup(dat_processed)

  return(dat_processed)
}
