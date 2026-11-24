#######################################################################
# Figure S7 - EBOV IgG optical density by serological and Plasmodium status (data)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
source("code/save_figure.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); d <- P$d
d$ebola_DO <- suppressWarnings(as.numeric(d$ebola_DO))
z <- d[complete.cases(d[, c("ebola_DO","ebola_IgG","paMalaria")]) & d$ebola_DO > 0, ]
draw <- function() {
par(mar = c(4, 4.5, 1, 1))
boxplot(log10(ebola_DO) ~ paMalaria + ebola_IgG, z,
        names = c("IgG- / Plasmodium-","IgG- / Plasmodium+","IgG+ / Plasmodium-","IgG+ / Plasmodium+"),
        ylab = "log10 EBOV IgG optical density", xlab = "", col = c("grey85","#f4a582","grey60","#d6604d"),
        outline = FALSE, las = 1, cex.axis = 0.7)
}
save_fig("FigS7_EBOV_optical_density", draw, width = 5.5, height = 4)
