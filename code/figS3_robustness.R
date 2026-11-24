#######################################################################
# Figure S3 - Key associations across model specifications (resultats/T11)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(ggplot2)}); source("code/save_figure.R", local = TRUE)
B <- read.csv("resultats/T11_key_associations_robustness.csv")
B$pair <- gsub("Ebola (IgG)", "EBOV IgG", paste(B$outcome, "/", B$exposure), fixed = TRUE)
long <- rbind(
  data.frame(pair = B$pair, spec = "Penalised GLMM, spline age (both directions)", OR = B$OR, lo = B$lo, hi = B$hi),
  data.frame(pair = B$pair, spec = "+ hunting, bushmeat, Pygmy community", OR = B$OR_behaviour, lo = NA, hi = NA),
  data.frame(pair = B$pair, spec = "+ all other pathogens", OR = B$OR_conditional, lo = NA, hi = NA),
  data.frame(pair = B$pair, spec = "Lagoon villages excluded", OR = B$OR_no_lagoon, lo = NA, hi = NA),
  data.frame(pair = B$pair, spec = "Leave-one-province-out (range)", OR = sqrt(B$LOPO_OR_min * B$LOPO_OR_max),
             lo = B$LOPO_OR_min, hi = B$LOPO_OR_max))
long$spec <- factor(long$spec, unique(long$spec))
long$pair <- factor(long$pair, rev(B$pair))
gA <- ggplot(long, aes(OR, pair, colour = spec)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", position = position_dodge(width = 0.7), na.rm = TRUE) +
  geom_point(position = position_dodge(width = 0.7), size = 1.6) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20)) +
  scale_colour_manual(values = c("black", "#e66101", "#5e3c99", "#1b9e77", "#b2abd2"), name = NULL) +
  geom_text(data = B, aes(x = 35, y = factor(pair, levels(long$pair)), label = sprintf("E = %.1f", Evalue_CI)),
            inherit.aes = FALSE, size = 2.4, hjust = 0) +
  coord_cartesian(xlim = c(0.7, 60)) +
  labs(x = "Odds ratio (log scale)", y = NULL) + theme_classic(base_size = 8) +
  theme(legend.position = "bottom", legend.direction = "vertical", axis.text.y = element_text(face = "italic"))
save_fig("FigS3_robustness", function() print(gA), width = 6.4, height = 5.4)
