#######################################################################
# 06 - Village level: village-level prevalence correlations and
#      habitat differences in pathogen community composition
# Manuscript items: Appendix S1 Section S8; Figure S6 (figS6_village_EBOV_Plasmodium.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab
vars <- c(P$patho12, "paMalaria")
nmin <- 15

vp <- do.call(rbind, lapply(split(d, d$village), function(z) {
  data.frame(village = z$village[1], ecosys2 = z$ecosys2[1], n = nrow(z),
             t(sapply(vars, function(v) if (sum(!is.na(z[[v]])) >= nmin) mean(z[[v]], na.rm = TRUE) else NA)))
}))
write.csv(vp, "resultats/T7_village_prevalence.csv", row.names = FALSE)

pairs <- list(c("ebola_IgG","paMalaria"), c("ebola_IgG","Pfalci"), c("Loaloa","paHTLV"),
              c("Loaloa","Fil_perstans"), c("dengue","westnile"), c("Pfalci","Pmalariae"),
              c("paMalaria","chik"), c("paMalaria","dengue"), c("paMalaria","paHTLV"))
vc <- do.call(rbind, lapply(pairs, function(pr) {
  z <- vp[complete.cases(vp[, pr]), ]
  ct <- suppressWarnings(cor.test(z[[pr[1]]], z[[pr[2]]], method = "spearman"))
  # partial correlation controlling for ecosystem (ranks residualised on ecosystem)
  r1 <- resid(lm(rank(z[[pr[1]]]) ~ z$ecosys2)); r2 <- resid(lm(rank(z[[pr[2]]]) ~ z$ecosys2))
  pc <- cor.test(r1, r2)
  data.frame(A = lab[pr[1]], B = lab[pr[2]], n_villages = nrow(z), rho = unname(ct$estimate),
             p = ct$p.value, rho_partial_ecosys = unname(pc$estimate), p_partial = pc$p.value)
}))
write.csv(vc, "resultats/T7_village_correlations.csv", row.names = FALSE)
print(vc, digits = 3, row.names = FALSE)

