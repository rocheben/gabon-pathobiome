#######################################################################
# Figure 4 - Association screening under three null models (resultats/T3_ASA_all.csv)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); lab <- P$lab
res <- read.csv("resultats/T3_ASA_all.csv")
fs <- res[res$set_id != "S8_all10", ]
fs$setlab <- sapply(strsplit(fs$set, "\\+"), function(v) paste(lab[v], collapse = " / "))
fs$null <- factor(fs$null, c("global","village","vill_age"),
                  c("Global independence","Within village","Within village x age"))
g <- ggplot(fs, aes(sub("^c", "", combo), obs / pmax(expected, 0.5), colour = null, shape = signif_envelope)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(position = position_dodge(width = 0.6), size = 1.3) +
  scale_y_log10() + scale_shape_manual(values = c(1, 16), labels = c("No", "Yes"), name = "Outside 95% envelope") +
  scale_colour_manual(values = c("grey40", "#1f78b4", "#e31a1c"), name = "Null model") +
  facet_wrap(~setlab, scales = "free_x", ncol = 2) +
  labs(x = "Combination (1 = positive, order as in panel title)", y = "Observed / expected (log)") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5),
                                  legend.position = "inside", legend.position.inside = c(0.76, 0.1),
                                  legend.box = "vertical", legend.box.just = "left", strip.text = element_text(face = "italic"))
save_fig("Fig4_association_screening", function() print(g), width = 6.4, height = 8.4)
