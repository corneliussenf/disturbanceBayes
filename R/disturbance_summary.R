#' Function for summarizing disturbances from a TimeSync sample
#'
#' @param dat Data exported from TimeSync.
#' @param by.agent Logical. Whether disturbances should be aggregated by agent.
#' @param grouping.vars Data.frame. Additional grouping variables stored in a data frame, where each row corresponds to one plot. The data frame must have a column 'plotid'.
#' @param change.agents Change agents to be considered. Standard is 'Harvest', 'Wind', and 'Decline'.
#' @return A sumary of the disturbances per year (and agents/grouping variables).
#' @export

disturbance_summary <- function(dat,
                                by.year = TRUE,
                                by.agent = FALSE,
                                grouping.vars = NULL,
                                sub.agents = FALSE,
                                change.agents = c("Harvest", "Wind", "Decline", "Hydrology", "Debris", "Other"),
                                agent.regroup = NULL,
                                interpreter = NULL,
                                remove.start = TRUE) {

  if (!is.null(interpreter)) dat <- dat[dat$interpreter == interpreter, ]

  if (sub.agents) {
    dat$change_process <- as.character(dat$change_process)
    agent_new <- unlist(lapply(strsplit(as.character(dat$change_process_notes), "|", fixed = TRUE), function(x) x[1]))
    dat[!is.na(agent_new) & dat$change_process == "Harvest", "change_process"] <- agent_new[!is.na(agent_new) & dat$change_process == "Harvest"]

    still_harvest <- dat[dat$change_process == "Harvest", "plotid"]
    print(paste0("Warning: Detected 'harvests' without sub-agent, which were automatically set to 'clearcut': ", paste(still_harvest, collapse = ", ")))
    dat[dat$change_process == "Harvest", "change_process"] <- "clearcut"

    change.agents <- c(change.agents, unique(agent_new))
    change.agents <- change.agents[!is.na(change.agents)]
  }

  dat_processed <- dplyr::mutate(dat, agent = lead(change_process))
  dat_processed <- dplyr::filter(dat_processed, agent %in% change.agents & dominant_landuse == "Forest")
  #dat_processed <- dplyr::select(dat_processed, -change_process)

  dat_processed$agent <- as.factor(as.character(dat_processed$agent))

  dat_processed <- dplyr::left_join(dat_processed, agent.regroup, by = "agent")
  dat_processed <- dplyr::mutate(dat_processed, "agent" = new)
  dat_processed <- dplyr::select(dat_processed, -new)

  if (is.null(grouping.vars)) {
    dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, "image_year", "agent"), disturbance = length(agent))
    forest <- sum(dplyr::summarize(dplyr::group_by(dat, plotid),
                                   forest = sum(dominant_landuse == "Forest") > 0)$forest)
    dat_processed <- dplyr::mutate(dat_processed, forest = forest)
  } else {

    dat_processed <- dplyr::left_join(dat_processed, grouping.vars, by = "plotid")
    grouping.vars.names <- names(grouping.vars)[-which(names(grouping.vars) == "plotid")]

    dat <- dplyr::left_join(dat, grouping.vars, by = "plotid")
    forest <- dplyr::summarize(dplyr::group_by_(dat, .dots = c("plotid", lapply(grouping.vars.names, function(x) x ))),
                                   forest = sum(dominant_landuse == "Forest") > 0)
    forest <- dplyr::summarize(dplyr::group_by_(forest, .dots = c(lapply(grouping.vars.names, function(x) x ))),
                               forest = sum(forest))

    dat_processed <- dplyr::summarise(dplyr::group_by_(dat_processed, .dots = c("image_year", lapply(grouping.vars.names, function(x) x ), "agent")), disturbance = length(agent))
    dat_processed <- dplyr::left_join(dat_processed, forest, by = grouping.vars.names)
  }

  dat_processed <- dplyr::mutate(dat_processed, year = image_year + 1)
  dat_processed <- dplyr::ungroup(dat_processed)
  dat_processed <- dplyr::select(dat_processed, -image_year)
  dat_processed$agent <- tolower(dat_processed$agent)
  #dat_processed <- dplyr::filter(dat_processed, !(agent == "decline" & year == 1985))
  if (remove.start == TRUE) dat_processed <- dplyr::filter(dat_processed, year > 1985)

  if (by.agent == FALSE) {
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

  dat_processed <- ungroup(dat_processed)

  return(dat_processed)
}
