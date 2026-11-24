#######################################################################
# Figure 1 - Structure of the pathogen community
#  (a) prevalence by ecosystem (resultats/T_prevalence_by_ecosystem.csv)
#  (b) number of pathogen groups per host vs independence (resultats/T_richness_distribution.csv)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); lab <- P$lab
pe <- read.csv("resultats/T_prevalence_by_ecosystem.csv")
pe$lo <- mapply(function(p, n) binom.test(round(p * n), n)$conf.int[1], pe$prev, pe$n)
pe$hi <- mapply(function(p, n) binom.test(round(p * n), n)$conf.int[2], pe$prev, pe$n)
pe$label <- factor(pe$label, lab[P$patho10])
pe$ecosys <- factor(pe$ecosys, c("FORET","SAVANE","LAGUNE"), c("Forest","Savanna","Lagoon"))
gA <- ggplot(pe, aes(label, prev, colour = ecosys)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5), size = 0.3) +
  scale_colour_manual(values = c("#1b7837", "#d9a520", "#2166ac"), name = NULL) +
  labs(x = NULL, y = "Prevalence / seroprevalence", tag = "a") + theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "top")
richness <- read.csv("resultats/T_richness_distribution.csv")
gB <- ggplot(richness, aes(k)) +
  geom_col(aes(y = observed), fill = "grey70") +
  geom_line(aes(y = expected_indep, colour = "Independence (global prevalence)"), linewidth = 0.8) +
  geom_point(aes(y = expected_indep, colour = "Independence (global prevalence)")) +
  geom_pointrange(aes(y = expected_village, ymin = village_lo, ymax = village_hi,
                      colour = "Independence within village"), position = position_nudge(x = 0.15)) +
  scale_colour_manual(values = c("#1f5f8b", "#c0392b"), name = NULL) +
  labs(x = "Number of pathogen groups per host", y = "Hosts", tag = "b") +
  scale_x_continuous(breaks = 0:10) + theme_classic(base_size = 9) + theme(legend.position = "top")
save_fig("Fig1_community_structure", function() {
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(1.1, 1))))
  print(gA, vp = grid::viewport(layout.pos.row = 1)); print(gB, vp = grid::viewport(layout.pos.row = 2))
}, width = 6.4, height = 6.4)
