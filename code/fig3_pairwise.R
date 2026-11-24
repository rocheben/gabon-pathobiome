#######################################################################
# Figure 3 - Pairwise associations, crude and adjusted (resultats/T2_pairwise.csv)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); lab <- P$lab
res <- read.csv("resultats/T2_pairwise.csv")
v12 <- P$patho12
hm <- rbind(
  data.frame(x = res$A, y = res$B, lor = log(res$OR_crude), q = res$q_crude, panel = "Crude"),
  data.frame(x = res$A, y = res$B, lor = log(res$OR_adj),  q = res$q_adj,  panel = "Adjusted (age, sex, village)"))
hm <- hm[hm$x %in% v12 & hm$y %in% v12, ]
hm$x <- factor(lab[hm$x], lab[v12]); hm$y <- factor(lab[hm$y], lab[v12])
hm$lor <- pmax(pmin(hm$lor, 2), -2)
hm$star <- ifelse(hm$q < 0.001, "***", ifelse(hm$q < 0.01, "**", ifelse(hm$q < 0.05, "*", "")))
hm$panel <- factor(hm$panel, c("Crude", "Adjusted (age, sex, village)"))
g <- ggplot(hm, aes(x, y, fill = lor)) + geom_tile(colour = "white") +
  geom_text(aes(label = star), size = 2.4) + facet_wrap(~panel) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-2, 2), name = "log OR", guide = guide_colourbar(raster = FALSE)) +
  labs(x = NULL, y = NULL) + theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), axis.text.y = element_text(face = "italic"))
save_fig("Fig3_pairwise", function() print(g), width = 6.4, height = 3.4)
