#' Function for estimating disturbance rates from TimeSync samples using Bayesian hierarchical modeling
#'
#' @param formula A model formula.
#' @param data A data.frame containing the data.
#' @param family The model family. Currently implemented are 'binomial' and 'poisson'.
#' @param summary_probs Vector containing the quantiles for summarizing the joint posterior distribution.
#' @return A sumary of the joint posterior distribution.
#' @export

bayes_estimator <- function(formula,
                            data,
                            family,
                            summary_probs = c(0.1, 0.5, 0.9)) {

  # Test
  #formula <- cbind(disturbance, forest - disturbance) ~ (1 | year : stratum)
  #formula <- cbind(disturbance, forest - disturbance) ~ (1 | year : stratum)
  #formula <- disturbance ~ (1 | year : stratum)
  #formula <- cbind(disturbance, forest - disturbance) ~ (1 | year) + (1 | stratum)
  #formula <- cbind(disturbance, forest - disturbance) ~ (1 | stratum / year)
  #formula <- cbind(disturbance, forest - disturbance) ~ (1 | year / stratum)
  #data <- sample_dat2
  #family <- "binomial"

  invlogit <- plogis  # function(x) 1/(1 + exp(-x))

  shift_draws <- function(draws) {
    sweep(draws[, -1], MARGIN = 1, STATS = draws[, 1], FUN = "+")
  }

  summary_stats <- function(posterior, f = family) {
    if (f == "binomial") {
      x <- invlogit(posterior)
    } else if (f == "poisson") {
      x <- posterior
    }

    t(apply(x, 2, quantile, probs = summary_probs))
  }

  if (family == "binomial") {

    ### Binomial ###

    vars.response <- all.vars(formula)[1:2]
    vars.predictor <- all.vars(formula)[3:length(all.vars(formula))]

    fit_partialpool <- rstanarm::stan_glmer(formula,
                                            data = data,
                                            family = stats::binomial("logit"))

    alphas <- shift_draws(as.matrix(fit_partialpool))
    partialpool <- summary_stats(alphas)
    partialpool <- as.data.frame(partialpool[-nrow(partialpool),])
    colnames(partialpool) <- c("lower", "estimate", "upper")

    if (length(vars.predictor) == 1) {

      random1 <- substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))
      random1_name <- unique(unlist(lapply(strsplit(random1, ":"), function(x) x[1])))
      partialpool[,random1_name] <- unlist(lapply(strsplit(random1, ":"), function(x) x[2]))
      rownames(partialpool) <- NULL

    } else if (length(vars.predictor) == 2) {

      if (length(attributes(stats::terms(formula))$term.labels) == 2) {

        partialpool <- partialpool[-nrow(partialpool), ]
        partialpool1 <- partialpool[grep(vars.predictor[1], substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))), ]
        partialpool2 <- partialpool[grep(vars.predictor[2], substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))), ]

        random1 <- substring(rownames(partialpool1), 15, (nchar(rownames(partialpool1)) - 1))
        random1_name <- unique(unlist(lapply(strsplit(random1, ":"), function(x) x[1])))
        partialpool1[,random1_name] <- unlist(lapply(strsplit(random1, ":"), function(x) x[2]))
        rownames(partialpool1) <- NULL

        random2 <- substring(rownames(partialpool2), 15, (nchar(rownames(partialpool2)) - 1))
        random2_name <- unique(unlist(lapply(strsplit(random2, ":"), function(x) x[1])))
        partialpool2[,random2_name] <- unlist(lapply(strsplit(random2, ":"), function(x) x[2]))
        rownames(partialpool2) <- NULL

        partialpool <- list()
        partialpool[[random1_name]] <- partialpool1
        partialpool[[random2_name]] <- partialpool2

      } else if (grepl("/", attributes(terms(formula))$term.labels, fixed = TRUE)) {

        partialpool <- partialpool[-nrow(partialpool), ]
        partialpool1 <- partialpool[grep(vars.predictor[2], substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))), ]
        partialpool2 <- partialpool[-grep(vars.predictor[2], substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))), ]

        random1a <- substring(rownames(partialpool1), 15, (nchar(rownames(partialpool1)) - 1))
        random1a_name <- unique(unlist(lapply(strsplit(random1a, ":"), function(x) x[1])))
        partialpool1[,random1a_name] <- unlist(lapply(strsplit(random1a, ":"), function(x) x[3]))

        random2a <- substring(rownames(partialpool1), 15, (nchar(rownames(partialpool1)) - 1))
        random2a_name <- unique(unlist(lapply(strsplit(random2a, ":"), function(x) x[2])))
        partialpool1[,random2a_name] <- unlist(lapply(strsplit(random2a, ":"), function(x) x[4]))

        rownames(partialpool1) <- NULL

        random2 <- substring(rownames(partialpool2), 15, (nchar(rownames(partialpool2)) - 1))
        random2_name <- unique(unlist(lapply(strsplit(random2, ":"), function(x) x[1])))
        partialpool2[,random2_name] <- unlist(lapply(strsplit(random2, ":"), function(x) x[2]))
        rownames(partialpool2) <- NULL

        partialpool <- list()
        partialpool[[paste(vars.predictor, collapse = ":")]] <- partialpool1
        partialpool[[random2_name]] <- partialpool2

      } else {

        random1 <- substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))
        random1_name <- unique(unlist(lapply(strsplit(random1, ":"), function(x) x[1])))
        partialpool[,random1_name] <- unlist(lapply(strsplit(random1, ":"), function(x) x[3]))

        random2 <- substring(rownames(partialpool), 15, (nchar(rownames(partialpool)) - 1))
        random2_name <- unique(unlist(lapply(strsplit(random2, ":"), function(x) x[2])))
        partialpool[,random2_name] <- unlist(lapply(strsplit(random2, ":"), function(x) x[4]))

        rownames(partialpool) <- NULL
      }

    }

    partialpool_mean <- as.data.frame(t(quantile(invlogit(as.matrix(fit_partialpool)[,1]), summary_probs)))
    colnames(partialpool_mean) <- c("lower", "estimate", "upper")

  } else if (family == "poisson") {

    ### Poisson ###

    vars.response <- all.vars(formula)[1]
    vars.predictor <- all.vars(formula)[2:length(all.vars(formula))]

    fit_partialpool <- rstanarm::stan_glmer(formula,
                                            data = data,
                                            family = stats::poisson("log"))

    alphas <- shift_draws(as.matrix(fit_partialpool))
    alphas <- alphas[,-ncol(alphas)]
    alphas <- data.table::melt(as.data.frame(alphas))

    if (length(vars.predictor) == 1) {
      random1 <- substring(as.character(alphas$variable), 15, (nchar(as.character(alphas$variable)) - 1))
      random1_name <- unique(unlist(lapply(strsplit(random1, ":"), function(x) x[1])))
      alphas[,random1_name] <- unlist(lapply(strsplit(random1, ":"), function(x) x[2]))
    } else if (length(vars.predictor) == 2) {
      random1 <- substring(as.character(alphas$variable), 15, (nchar(as.character(alphas$variable)) - 1))
      random1_name <- unique(unlist(lapply(strsplit(random1, ":"), function(x) x[1])))
      alphas[,random1_name] <- unlist(lapply(strsplit(random1, ":"), function(x) x[3]))

      random2 <- substring(as.character(alphas$variable), 15, (nchar(as.character(alphas$variable)) - 1))
      random2_name <- unique(unlist(lapply(strsplit(random2, ":"), function(x) x[2])))
      alphas[,random2_name] <- unlist(lapply(strsplit(random2, ":"), function(x) x[4]))
    }

    alphas$variable <- gsub("[^0-9]", "", alphas$variable)
    alphas <- dplyr::left_join(alphas, data, by = vars.predictor)
    alphas$dr <- exp(alphas$value) / alphas$forest

    if (length(vars.predictor) == 1) {

      partialpool <- dplyr::summarize(dplyr::group_by_(alphas, vars.predictor[1]),
                               lower = stats::quantile(dr, summary_probs[1]),
                               estimate = stats::quantile(dr, summary_probs[2]),
                               upper = stats::quantile(dr, summary_probs[3]))

    } else if (length(vars.predictor) == 2) {

      partialpool <- dplyr::summarize(dplyr::group_by_(alphas, vars.predictor[1], vars.predictor[2]),
                               lower = stats::quantile(dr, summary_probs[1]),
                               estimate = stats::quantile(dr, summary_probs[2]),
                               upper = stats::quantile(dr, summary_probs[3]))

    }

    partialpool_mean <- shift_draws(as.matrix(fit_partialpool))
    partialpool_mean <- exp(partialpool_mean[,ncol(partialpool_mean)]) / sample_dat[1, "forest"]
    partialpool_mean <- as.data.frame(t(quantile(partialpool_mean, summary_probs)))

  }

  return(list(random = partialpool, fixed = partialpool_mean, model = fit_partialpool))

}
