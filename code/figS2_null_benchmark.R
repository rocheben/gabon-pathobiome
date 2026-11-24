#######################################################################
# Figure S2 - Benchmark of null models: (a) type I error (resultats/T12), (b) power (resultats/T13)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
S <- read.csv("resultats/T12_ASA_null_simulation.csv"); W <- read.csv("resultats/T13_ASA_null_power.csv")
lv <- c("global","village","vill_age"); ll <- c("Global independence","Within village","Within village x age class")
S$null <- factor(S$null, lv, ll); W$null <- factor(W$null, lv, ll)
S$set <- factor(S$set, c("S1","S5"), c("L. loa / M. perstans /\nPlasmodium", "Plasmodium / EBOV IgG /\nCHIKV"))
cols <- c("grey35", "#1f78b4", "#e31a1c")
gA <- ggplot(S, aes(set, false_positive_rate_envelope, fill = null)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) + geom_hline(yintercept = 0.05, linetype = 2) +
  scale_fill_manual(values = cols, name = "Null model") +
  labs(x = "Pathogen set (simulated independence)", y = "False-positive rate", tag = "a") +
  theme_classic(base_size = 8) + theme(legend.position = "none", axis.text.x = element_text(face = "italic"))
gB <- ggplot(W, aes(OR, power, colour = null)) + geom_line() + geom_point(size = 1.8) +
  scale_colour_manual(values = cols, name = "Null model") + scale_x_continuous(breaks = c(1.5, 2, 3)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Injected odds ratio (Plasmodium / EBOV IgG)", y = "Power", tag = "b") +
  theme_classic(base_size = 8) + theme(legend.position = "inside", legend.position.inside = c(0.68, 0.22))
save_fig("FigS2_null_benchmark", function() {
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(gA, vp = grid::viewport(layout.pos.col = 1)); print(gB, vp = grid::viewport(layout.pos.col = 2))
}, width = 6.4, height = 3.4)
