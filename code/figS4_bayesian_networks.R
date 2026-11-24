#######################################################################
# Figure S4 - Model-averaged Bayesian networks, 10 pathogen groups, without and
#             with environmental and host nodes (resultats/T6_BN_arc_strength.csv)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(igraph)); source("code/save_figure.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); lab <- P$lab; p10 <- P$patho10
tab <- read.csv("resultats/T6_BN_arc_strength.csv")
draw <- function() {
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1), cex.main = 0.9)
lay <- NULL
for (s in c("10 taxa", "10 taxa + environment/host")) {
  t <- tab[tab$set == s & tab$in_avg, ]
  g <- graph_from_data_frame(t[, c("from","to")], directed = FALSE, vertices = data.frame(name = p10))
  if (is.null(lay)) lay <- layout_in_circle(g)
  E(g)$color <- ifelse(t$sign > 1, "#b2182b", "#2166ac")
  E(g)$width <- 1 + 6 * (t$strength - 0.5)
  gu <- c(Loaloa = "#e6ab02", Fil_perstans = "#e6ab02", paMalaria = "#66a61e")
  V(g)$color <- ifelse(V(g)$name %in% names(gu), gu[V(g)$name], "#7fa7d9")
  plot(g, layout = lay, vertex.label = lab[V(g)$name], vertex.label.cex = 0.55,
       vertex.label.color = "black", vertex.size = 24, vertex.label.font = 3,
       main = ifelse(s == "10 taxa", "Pathogens only", "With province, ecosystem,\nage and sex"))
}
}
save_fig("FigS4_bayesian_networks", draw, width = 6.4, height = 3.4)
