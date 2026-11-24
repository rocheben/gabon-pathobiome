#######################################################################
# Figure 2 - Scales of structure: (a) Hmsc variance partitioning (resultats/T14),
#            (b) village- vs individual-scale associations (resultats/T15)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
long <- read.csv("resultats/T14_hmsc_varpart_summary.csv")
long$component <- factor(long$component, c("Ecosystem", "Village (latent)", "Host (age, sex)", "Individual (latent)"))
ord <- long[long$component == "Individual (latent)", ]; ord <- ord$pathogen[order(ord$share)]
long$pathogen <- factor(long$pathogen, ord)
gA <- ggplot(long, aes(share, pathogen, fill = component)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c("#1b7837", "#a6dba0", "#bababa", "#762a83"), name = NULL) +
  labs(x = "Share of explained variation", y = NULL, tag = "a") + theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "right") +
  guides(fill = guide_legend(ncol = 1))
sc <- read.csv("resultats/T15_scale_concordance.csv")
sc$lab <- gsub("Ebola (IgG)", "EBOV IgG", ifelse(sc$class != "Neither", paste(sc$A, sc$B, sep = " / "), ""), fixed = TRUE)
sc$vj <- ifelse(grepl("L. loa / HTLV", sc$lab), 1.8, -0.8)
gB <- ggplot(sc, aes(rho_village, logOR_individual, colour = class)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 0, linetype = 3) +
  geom_point(size = 1.8) +
  geom_text(aes(label = lab, vjust = vj), size = 2.6, show.legend = FALSE, fontface = "italic") +
  scale_colour_manual(values = c(`Both scales` = "#762a83", `Village scale only` = "#1b7837",
                                 `Individual scale only` = "#e66101", Neither = "grey70"), name = NULL) +
  labs(x = "Village-scale association (partial Spearman rho)", y = "Individual-scale association (adjusted log OR)", tag = "b") +
  scale_x_continuous(limits = c(-0.35, 0.45)) + theme_classic(base_size = 9) + theme(legend.position = "right")
save_fig("Fig2_scales", function() {
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(1, 1.15))))
  print(gA, vp = grid::viewport(layout.pos.row = 1)); print(gB, vp = grid::viewport(layout.pos.row = 2))
}, width = 6.4, height = 7.2)
