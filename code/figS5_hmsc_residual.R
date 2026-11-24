#######################################################################
# Figure S5 - Individual-level residual associations, Hmsc model with nine pathogen
#             groups (resultats/T9_hmsc_associations_9.csv); only supported associations coloured
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
source("code/save_figure.R", local = TRUE)
h <- read.csv("resultats/T9_hmsc_associations_9.csv"); h <- h[h$level == "individual", ]
sp <- unique(c(h$A, h$B)); k <- length(sp)
M <- matrix(NA, k, k, dimnames = list(sp, sp))
for (i in seq_len(nrow(h))) { v <- ifelse(h$supported[i], h$cor_mean[i], 0); M[h$A[i], h$B[i]] <- v; M[h$B[i], h$A[i]] <- v }
draw <- function() {
par(mar = c(6.5, 6.5, 1, 1))
image(1:k, 1:k, M[, k:1], zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
      col = colorRampPalette(c("#2166ac", "white", "#b2182b"))(41))
axis(1, 1:k, sp, las = 2, cex.axis = 0.6, font = 3); axis(2, 1:k, rev(sp), las = 2, cex.axis = 0.6, font = 3)
}
save_fig("FigS5_hmsc_residual", draw, width = 4.5, height = 4.2)
