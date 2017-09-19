
library(raster)
library(landscapeR)
library(rstanarm)
library(tidyverse)

#devtools::install_github("corneliussenf/disturbanceBayes")
#library(disturbanceBayes)

source("R/bayes_estimator_unpooled.R")
source("R/bayes_estimator.R")
source("R/disturbance_summary.R")
#source("../mapping/lib/misc/grid_plots_shared_legend.R")

# Create random landscape -------------------------------------------------

ras <- raster(ncol = 1000, nrow = 1000, xmn = 1, xmx = 1000, ymn = 1, ymx = 1000)
values(ras) <- 0

landscape <- makeClass(ras, npatch = 250, size = round(runif(250, 5, 3500), 0))
fcr <- cellStats(landscape == 1, mean, na.rm = TRUE)
i <- 2
while(i < 32) {
  n <- floor(runif(1, 100, 1000))
  s <- floor(rlnorm(n, log(5), log(3)))
  s <- ifelse(s == 0, 1, s)
  tryCatch(landscape <- makeClass(landscape, npatch = n, size = s, val = i, bgr = 1), error = function(e) print("Something fishy going on..."))
  i <- i + 1
}

landscape[landscape == 0] <- NA

dist_rate_true <- table(as.vector(landscape))
dist_rate_true <- data.frame(disturbance = as.integer(dist_rate_true[2:length(dist_rate_true)]),
                             forest = as.integer(dist_rate_true[1]))

dist_rate_true$disturbance_rate <- dist_rate_true$disturbance / dist_rate_true$forest
dist_rate_true$year <- factor(1:30)

dist_rate_true_mean <- mean(dist_rate_true$disturbance_rate)

# Add agent

patches <- (landscape > 1) %>%
  clump(.)
patchid <- unique(na.omit(values(patches)))
agents <- data.frame(as.integer(patchid),
                     sample(c(1, 2, 3, 4), length(patchid), replace = TRUE, prob = c(0.5, 0.3, 0.2, 0.2)))
agents <- reclassify(patches, as.matrix(agents))
agents[is.na(agents)] <- 0
agents <- mask(agents, landscape)

landscape_agent <- na.omit(as.data.frame(stack(landscape, agents)))
landscape_agent_disturbance <- landscape_agent[landscape_agent[,1] > 1, ]

dist_rate_true_agent <- data.table::melt(table(landscape_agent_disturbance[,1], landscape_agent_disturbance[,2]))
names(dist_rate_true_agent) <- c("year", "agent", "disturbance")
dist_rate_true_agent$year <- factor(dist_rate_true_agent$year - 1)
dist_rate_true_agent$agent <- factor(dist_rate_true_agent$agent)
dist_rate_true_agent$forest <- sum(landscape_agent[,1] == 1)
dist_rate_true_agent$disturbance_rate <- dist_rate_true_agent$disturbance / dist_rate_true_agent$forest

# Add stratum

xy <- as.data.frame(ras, xy = TRUE)
xy_clust <- kmeans(xy[,1:2], centers = 25)
xy$stratum <- xy_clust$cluster
stratum <- rasterFromXYZ(xy[,-3])

landscape_stratum <- na.omit(as.data.frame(stack(landscape, stratum)))
landscape_stratum_disturbance <- landscape_stratum[landscape_stratum[,1] > 1, ]

dist_rate_true_stratum <- data.table::melt(table(landscape_stratum_disturbance[,1], landscape_stratum_disturbance[,2]))
names(dist_rate_true_stratum) <- c("year", "stratum", "disturbance")
dist_rate_true_stratum$year <- factor(dist_rate_true_stratum$year - 1)
dist_rate_true_stratum <- landscape_stratum %>%
  group_by(stratum) %>%
  summarize(forest = sum(layer == 1)) %>%
  right_join(dist_rate_true_stratum, by = "stratum")
dist_rate_true_stratum$stratum <- factor(dist_rate_true_stratum$stratum)
dist_rate_true_stratum$disturbance_rate <- dist_rate_true_stratum$disturbance / dist_rate_true_stratum$forest

dist_rate_true_stratum_2 <- dist_rate_true_stratum %>%
  group_by(stratum) %>%
  dplyr::summarize(disturbance_rate = mean(disturbance_rate))

# Agent and stratum

landscape_stratum_agent <- na.omit(as.data.frame(stack(landscape, stratum, agents)))
landscape_stratum_agent_disturbance <- landscape_stratum_agent[landscape_stratum_agent[,1] > 1, ]

dist_rate_true_stratum_agent <- data.table::melt(table(landscape_stratum_agent_disturbance[,2],
                                                       landscape_stratum_agent_disturbance[,3]))
names(dist_rate_true_stratum_agent) <- c("stratum", "agent", "disturbance")
dist_rate_true_stratum_agent <- landscape_stratum_agent %>%
  group_by(stratum) %>%
  summarize(forest = sum(layer.1 == 1)) %>%
  right_join(dist_rate_true_stratum_agent, by = "stratum")
dist_rate_true_stratum_agent$stratum <- factor(dist_rate_true_stratum_agent$stratum)
dist_rate_true_stratum_agent$agent <- factor(dist_rate_true_stratum_agent$agent)
dist_rate_true_stratum_agent$disturbance_rate <- dist_rate_true_stratum_agent$disturbance / dist_rate_true_stratum_agent$forest

# Plot simulated maps

library(RColorBrewer)
library(rasterVis)

getPalette <- colorRampPalette(brewer.pal(9, "BuPu"))

r <- ratify(landscape)
rat <- levels(r)[[1]]
rat$landcover <- c("Forest", 1985:2014)
levels(r) <- rat

p_year <- levelplot(r, col.regions = c(getPalette(30), brewer.pal(5, name = "Set1")[3]),
                    contour = FALSE, margin = FALSE, main = list("a)", x = 0.05),
                    xlab =  NULL, ylab =NULL, scales = list(draw = FALSE),
                    colorkey = list(at = (0:31) + 0.5,
                                    labels = list(labels = c("1985", rep("", 4),
                                               "1990", rep("", 4),
                                               "1995", rep("", 4),
                                               "2000", rep("", 4),
                                               "2005", rep("", 4),
                                               "2010", rep("", 4), "Forest"), cex = 0.75)))

r <- ratify(agents)
rat <- levels(r)[[1]]
rat$landcover <- c("Forest", "A", "B", "C", "D")
#rat$class <- c('A1', 'B2', 'C3')
levels(r) <- rat

p_agent <- levelplot(r, col.regions = RColorBrewer::brewer.pal(5, name = "Set1")[c(1, 2, 5, 4, 3)],
                     contour = FALSE, margin = FALSE, main = list("b)", x = 0.05),
                     xlab =  NULL, ylab =NULL, scales = list(draw = FALSE),
                     colorkey = list(labels = list(cex = 0.75)))


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

r <- ratify(stratum * (landscape != 1))
rat <- levels(r)[[1]]
rat$landcover <- factor(c("Forest", 1:25), levels = c(1:25, "Forest"))
levels(r) <- rat

p_stratum <- levelplot(r, col.regions = c(sample(col_vector, 25), brewer.pal(5, name = "Set1")[3]),
                     contour = FALSE, margin = FALSE, main = list("c)", x = 0.05),
                     xlab = NULL, ylab =NULL, scales = list(draw = FALSE),
                     colorkey = list(labels = list(cex = 0.75)))

p <- grid.arrange(p_year, p_agent, p_stratum, ncol = 2)
ggsave("simulated_data.pdf", p, path = "analysis/figures/", width = 8, height = 7.5)

# Save maps

save(landscape, agents, stratum, file = "analysis/results/simulation_datasets.RData")
load("analysis/results/simulation_datasets.RData")

# Simulation study analysis -----------------------------------------------

### Sample sizes

sn <- 4000

### Function for analyzing models

analyze_models <- function(modellist, true, name, sim_dat) {

  # Join estimate and true

  join <- modellist %>%
    map(~ left_join(.$estimate, true, by = names(true)[which(!(names(true) %in% c("disturbance", "forest", "disturbance_rate")))]))

  # Calculate accuracy measures

  acc <- join %>%
    map(~ data.frame(model = name,
                     mae = mean(abs((.$Q50 * 100) - (.$disturbance_rate * 100))),
                     mae_rel = mean(abs((.$Q50 * 100) - (.$disturbance_rate * 100))) / mean(.$disturbance_rate * 100),
                     mse = mean(((.$Q50 * 100) - (.$disturbance_rate * 100))^2),
                     mse_rel = mean(((.$Q50 * 100) - (.$disturbance_rate * 100))^2) / mean(.$disturbance_rate * 100),
                     rmse = sqrt(mean(((.$Q50 * 100) - (.$disturbance_rate * 100))^2)),
                     r = cor((.$Q50 * 100), (.$disturbance_rate * 100)))) %>%
    set_names(c("Partially pooled", "Un-pooled")) %>%
    bind_rows(.id = "type")

  # Plots

  sim_dat <- sim_dat %>%
    rename(sample_n = disturbance) %>%
    rename(forest_n = forest)

  plotdat <- join %>%
    set_names(c("Partially pooled", "Un-pooled")) %>%
    bind_rows(.id = "type") %>%
    left_join(sim_dat) %>%
    mutate(id = rep(1:nrow(sim_dat), 2)) %>%
    mutate(difference = (Q50 * 100) - (disturbance_rate * 100),
           difference_abs = abs(difference))

  max_lim <- max(c(plotdat$Q50 * 100, plotdat$disturbance_rate * 100)) * 1.01

  scatter <- ggplot(plotdat, aes(y = Q50 * 100, x = disturbance_rate * 100, shape = type, col = type)) +
    geom_point() +
    xlim(0, max_lim) +
    ylim(0, max_lim) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(x = "True", y = "Estimated", col = NULL, shape = NULL) +
    ggthemes::theme_base() +
    scale_color_brewer(palette = "Set1") +
    theme(legend.justification = c(1, 0),
          legend.position = c(1, 0),
          plot.background = element_blank(),
          legend.background = element_blank())

  # density <- ggplot(plotdat, aes(x = (Q50 * 100) - (disturbance_rate * 100), col = type)) +
  #   stat_ecdf() +
  #   geom_vline(xintercept = 0) +
  #   labs(x = "Estimated minus true disturbance rate", y = "Empirical CDF", col = NULL) +
  #   ggthemes::theme_base() +
  #   scale_color_brewer(palette = "Set1") +
  #   theme(legend.position = "bottom",
  #         plot.background = element_blank())

  density <- ggplot(plotdat, aes(x = difference, y = as.integer(reorder(id, difference)))) +
    geom_point(aes(col = type, shape = type), alpha = 0.75)  +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Estimation error", y = "Estimate (orderd by estimation error)", col = NULL, shape = NULL) +
    ggthemes::theme_base() +
    scale_color_brewer(palette = "Set1") +
    theme(legend.justification = c(1, 0),
          legend.position = c(1, 0),
          plot.background = element_blank(),
          legend.background = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  return(list(accuracy = acc,
              plotdat = plotdat,
              scatterplot = scatter,
              density = density))

}

### Run models

### Estimate disturbance rate per agent : year

sample <- sampleRandom(stack(landscape, agents), size = sn, na.rm = TRUE)
sample_disturbance <- sample[sample[,1] > 1, ]

simulation_dat <- data.table::melt(table(sample_disturbance[,1], sample_disturbance[,2]))
names(simulation_dat) <- c("year", "agent", "disturbance")
simulation_dat$year <- factor(simulation_dat$year - 1)
simulation_dat$agent <- factor(simulation_dat$agent)
simulation_dat$forest <- sum(sample[,1] == 1)

partialpooled <- bayes_estimator(formula = cbind(disturbance, forest - disturbance) ~ (1 | agent : year),
                                          data = simulation_dat,
                                          family = "binomial")

unpooled <- bayes_estimator_unpooled(formula = cbind(disturbance, forest - disturbance) ~ agent : year,
                                              data = simulation_dat,
                                              family = "binomial")

models_year_agent <- list(partialpooled, unpooled)
results_year_agent <- analyze_models(modellist = models_year_agent,
                                     true = dist_rate_true_agent,
                                     name = "by_year_agent",
                                     sim_dat = simulation_dat)

### Estimate disturbance rate per stratum : agent

sample <- sampleRandom(stack(landscape, agents, stratum), size = sn, na.rm = TRUE)
sample_disturbance <- sample[sample[,1] > 1, ]

simulation_dat <- data.table::melt(table(sample_disturbance[,2], sample_disturbance[,3]))
names(simulation_dat) <- c("agent", "stratum", "disturbance")
simulation_dat <- sample %>%
  as.data.frame(.) %>%
  group_by(stratum) %>%
  summarize(forest = sum(layer.1 == 1)) %>%
  right_join(simulation_dat, by = "stratum")
simulation_dat$stratum <- factor(simulation_dat$stratum)
simulation_dat$agent <- factor(simulation_dat$agent)

partialpooled <- bayes_estimator(formula = cbind(disturbance, forest - disturbance) ~ (1 | stratum : agent),
                                          data = simulation_dat,
                                          family = "binomial")

unpooled <- bayes_estimator_unpooled(formula = cbind(disturbance, forest - disturbance) ~ stratum : agent,
                                              data = simulation_dat,
                                              family = "binomial")

models_agent_stratum <- list(partialpooled, unpooled)
results_year_stratum <- analyze_models(modellist = models_agent_stratum,
                                       true = dist_rate_true_stratum_agent,
                                       name = "by_agent_stratum",
                                       sim_dat = simulation_dat)

### Estimate disturbance rate per stratum : agent : year

sample <- sampleRandom(stack(landscape, agents, stratum), size = sn, na.rm = TRUE)
sample_disturbance <- sample[sample[,1] > 1, ]

simulation_dat <- data.table::melt(table(sample_disturbance[,1], sample_disturbance[,2], sample_disturbance[,3]))
names(simulation_dat) <- c("year", "agent", "stratum", "disturbance")
simulation_dat$year <- factor(simulation_dat$year - 1)
simulation_dat <- sample %>%
  as.data.frame(.) %>%
  group_by(stratum) %>%
  summarize(forest = sum(layer.1 == 1)) %>%
  right_join(simulation_dat, by = "stratum")
simulation_dat$stratum <- factor(simulation_dat$stratum)
simulation_dat$agent <- factor(simulation_dat$agent)

partialpooled <- bayes_estimator(formula = cbind(disturbance, forest - disturbance) ~ (1 | stratum : agent : year),
                                 data = simulation_dat,
                                 family = "binomial")

unpooled <- bayes_estimator_unpooled(formula = cbind(disturbance, forest - disturbance) ~ stratum : agent : year,
                                     data = simulation_dat,
                                     family = "binomial")

models_year_agent_stratum <- list(partialpooled, unpooled)
results_year_agent_stratum <- analyze_models(models_year_agent_stratum, dist_rate_true_stratum_agent, "by_year_stratum_agent")

### Save everything

save(models_year, models_year_agent, models_year_stratum, models_year_agent_stratum, file = "analysis/results/models.RData")

### Combine plots and create tables

# Plot

plots <- grid_arrange_shared_legend(list(results_year$scatterplot +
                                           labs(title = "Year") +
                                           theme(plot.title = element_text(size = 11),
                                                 axis.title = element_text(size = 11),
                                                 axis.text = element_text(size = 11),
                                                 legend.text = element_text(size = 11)),
                                         results_year_agent$scatterplot +
                                           labs(title = "Year-Agent") +
                                           theme(plot.title = element_text(size = 11),
                                                 axis.title = element_text(size = 11),
                                                 axis.text = element_text(size = 11)),
                                         results_year_stratum$scatterplot +
                                           labs(title = "Year-Stratum") +
                                           theme(plot.title = element_text(size = 11),
                                                 axis.title = element_text(size = 11),
                                                 axis.text = element_text(size = 11)),
                                         results_year$density +
                                           theme(axis.title = element_text(size = 11),
                                                 axis.text = element_text(size = 11)) +
                                           labs(title = "",
                                                x = "Estimated - True"),
                                         results_year_agent$density +
                                           theme(axis.title = element_text(size = 11),
                                                 axis.text = element_text(size = 11)) +
                                           labs(title = "",
                                                x = "Estimated - True"),
                                         results_year_stratum$density +
                                           theme(axis.title = element_text(size = 11),
                                                 axis.text = element_text(size = 11)) +
                                           labs(title = "",
                                                x = "Estimated - True")),
                                         # results_year_agent_stratum$scatterplot +
                                         #   labs(title = "Year-Agent-Stratum") +
                                         #   theme(plot.title = element_text(size = 12),
                                         #         axis.title = element_text(size = 12),
                                         #         axis.text = element_text(size = 10)),
                                         # results_year_agent_stratum$density +
                                         #   theme(axis.title = element_text(size = 12),
                                         #         axis.text = element_text(size = 10)) +
                                         #   labs(title = "",
                                         #        x = "Estimated - True")),
                                    ncol = 3, nrow = 2, position = "bottom")

grid.draw(plots)

ggsave("simulation_results.pdf", plots, path = "analysis/figures/", width = 8, height = 6)

# Table

list(results_year,
     results_year_agent,
     results_year_stratum,
     results_year_agent_stratum) %>%
  map(~ .$accuracy) %>%
  map(~ data.table::dcast(., model ~ type, value.var = "mse_rel")) %>%
  bind_rows() %>%
  write.csv(., "analysis/results/model_performance_rmse.csv", row.names = FALSE)

# Interpreter error -------------------------------------------------------



# With TimeSync data ------------------------------------------------------

# Load data

dat <- read.csv("analysis/data/austria_vertex_07_07_17.csv")

group <- read.csv("analysis/data/Samples mit WGB_BFI_und_Hoehe.csv")
group$plotid <- group$id

wuchsgebiete <- group[, c("plotid", "Wuchsge1")]
wuchsgebiete$Wuchsge1 <- as.character(wuchsgebiete$Wuchsge1)

hoehenstufen <- group[, c("plotid", "SRTMkomple")]
hoehenstufen$SRTMkomple <- cut(hoehenstufen$SRTMkomple,
                               seq(0, 4000, 400),
                               include.lowest = TRUE,
                               labels = paste0(seq(0, 3600, 400), "-", seq(400, 4000, 400)))

regroup <- data.frame(agent = c("Harvest", "Wind", "Decline", "Hydrology", "Debris", "Other"),
                      new = c("Human", rep("Natural", 5)))

wuchsgebiete_shape <- raster::shapefile("analysis/data/wuchsgebiete_shape/WLamPoly.shp")

# Grouped by Wuchsgebiet

disturbance_dat <- disturbance_summary(dat,
                                       by.year = FALSE,
                                       grouping.vars = wuchsgebiete,
                                       agent.regroup = regroup,
                                       by.agent = FALSE,
                                       interpreter = 131)

model1 <- bayes_estimator(formula = cbind(disturbance, forest - disturbance) ~ (1 | Wuchsge1),
                         data = disturbance_dat,
                         family = "binomial")

ggplot(model1$estimate,
       aes(x = reorder(Wuchsge1, Q50), y = Q50 * 100 / 32)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = Q025 * 100 / 32, ymax = Q975 * 100 / 32),
                position = position_dodge(width = 0.25), width = 0) +
  coord_flip() +
  labs(x = NULL, y = expression(paste("Disturbance rate [%yr"^"-1","]")), col = NULL, shape = NULL) +
  ggthemes::theme_base() +
  scale_color_manual(values = c("grey20", "grey60")) +
  theme(legend.justification = c(1, 0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        plot.background = element_blank())

wuchsgebiete_shape.fortify <- fortify(wuchsgebiete_shape, region = "Wuchsge1")
wuchsgebiete_shape.fortify <- left_join(wuchsgebiete_shape.fortify,
                                        model1$estimate,
                                        by = c("id" = "Wuchsge1"))
#wuchsgebiete_shape.fortify <- data.table::melt(wuchsgebiete_shape.fortify, measure.vars = c("human", "natural"))

p <- ggplot(wuchsgebiete_shape.fortify) +
  aes(long,lat,group = group, fill = Q50 * 100 / 32) +
  geom_polygon() +
  geom_path(color = "white") +
  coord_equal() +
  scale_fill_gradient(low = "black", high = "grey90") +
  labs(x = NULL, y = NULL, fill = expression(paste("Disturbance rate [%yr"^"-1","]"))) +
  ggthemes::theme_map() +
  theme(legend.justification = c(0, 1),
        legend.position = c(0, 1),
        legend.background = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5,
                               direction = "horizontal", title.position = "top"))

ggsave("austria_wuchsgebiete.pdf", p, path = "analysis/figures/", width = 3.5, height = 2)

# Grouped by year and agent

disturbance_dat <- disturbance_summary(dat,
                                       by.year = TRUE,
                                       agent.regroup = regroup,
                                       by.agent = FALSE,
                                       interpreter = 131)

#disturbance_dat$agent <- ifelse(disturbance_dat$agent == "human", "Human", "Natural")

model2 <- bayes_estimator(formula = cbind(disturbance, forest - disturbance) ~ (1 | year),
                         data = disturbance_dat,
                         family = "binomial")

p <- ggplot(model2$estimate,
       aes(x = as.integer(year), y = Q50 * 100)) +
  #geom_point(position = position_dodge(width = 1)) +
  #geom_errorbar(aes(ymin = Q025 * 100, ymax = Q975 * 100),
  #              position = position_dodge(width = 1), width = 0) +
  geom_ribbon(aes(ymin = Q025 * 100, ymax = Q975 * 100), alpha = 0.3, col = NA) +
  geom_line() +
  labs(x = "Year", y = expression(paste("Disturbance rate [%yr"^"-1","]")), col = NULL, shape = NULL) +
  ggthemes::theme_base() +
  scale_color_manual(values = c("grey20", "grey60")) +
  theme(legend.justification = c(0, 1),
        legend.position = c(0, 1),
        legend.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave("austria_year_agent.pdf", p, path = "analysis/figures/", width = 3.5, height = 3.5)

# Grouped by agent and elevation

disturbance_dat <- disturbance_summary(dat,
                                       by.year = FALSE,
                                       grouping.vars = hoehenstufen,
                                       agent.regroup = regroup,
                                       by.agent = FALSE,
                                       interpreter = 131)

#disturbance_dat$agent <- ifelse(disturbance_dat$agent == "human", "Human", "Natural")

model3 <- bayes_estimator(formula = cbind(disturbance, forest - disturbance) ~ (1 | SRTMkomple),
                          data = disturbance_dat,
                          family = "binomial")

model3$estimate$SRTMkomple <- factor(model3$estimate$SRTMkomple, levels = paste0(seq(0, 3600, 400), "-", seq(400, 4000, 400)))

p <- ggplot(model3$estimate,
       aes(x = SRTMkomple, y = Q50 * 100 / 32)) +
  geom_point(position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = Q025 * 100 / 32, ymax = Q975 * 100 / 32),
                position = position_dodge(width = 0.25), width = 0) +
  coord_flip() +
  labs(x = "Elevation [m]", y = expression(paste("Disturbance rate [%yr"^"-1","]")), col = NULL, shape = NULL) +
  ggthemes::theme_base() +
  scale_color_manual(values = c("grey20", "grey60")) +
  theme(legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.background = element_blank(),
        plot.background = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave("austria_elevation.pdf", p, path = "analysis/figures/", width = 3.5, height = 3.5)
