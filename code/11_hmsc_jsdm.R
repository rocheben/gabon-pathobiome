#######################################################################
# 11 - Joint species distribution model (Hmsc, Ovaskainen et al. 2017)
#  Probit JSDM of pathogen occurrence with age (quadratic), sex and
#  ecosystem as fixed effects, and latent factors at two levels:
#  village (shared exposure not captured by covariates) and individual
#  (residual within-host co-occurrence). Individual-level residual
#  associations with >= 95% posterior support are the JSDM analogue of
#  "associations not explained by shared exposure".
# Manuscript items: Tables S6 and S7; Figure S5 (figS5_hmsc_residual.R); Figure 2a (via 12_scales.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(Hmsc)})
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab
for (set in c("9", "10")) {
vars <- if (set == "9") setdiff(P$patho10, "westnile") else P$patho10
dd <- d[complete.cases(d[, c(vars, "age", "sex01", "ecosys2", "village")]), ]
Y <- as.matrix(dd[, vars])
XData <- data.frame(age = as.numeric(scale(dd$age)), sex = dd$sex01, eco = droplevels(dd$ecosys2))
studyDesign <- data.frame(village = factor(dd$village), indiv = factor(seq_len(nrow(dd))))
rlV <- HmscRandomLevel(units = studyDesign$village)
rlI <- HmscRandomLevel(units = studyDesign$indiv)
m <- Hmsc(Y = Y, XData = XData, XFormula = ~ age + I(age^2) + sex + eco, distr = "probit",
          studyDesign = studyDesign, ranLevels = list(village = rlV, indiv = rlI))
set.seed(11)
t0 <- Sys.time()
message("Hmsc MCMC, ", set, " pathogen groups (no progress bar; several minutes)...")
m <- sampleMcmc(m, samples = 1000, thin = 10, transient = 5000, nChains = 2, nParallel = 2, verbose = 0)
cat("MCMC time:", format(Sys.time() - t0), "\n")
saveRDS(m, paste0("resultats/hmsc_model_", set, ".rds"))

post <- convertToCodaObject(m)
psrf_beta <- gelman.diag(post$Beta, multivariate = FALSE)$psrf[, 1]
OmegaI <- computeAssociations(m)[[2]]; OmegaV <- computeAssociations(m)[[1]]
tab <- do.call(rbind, lapply(1:2, function(k) {
  A <- computeAssociations(m)[[k]]; nm <- c("village","individual")[k]
  idx <- which(upper.tri(A$mean), arr.ind = TRUE)
  data.frame(level = nm, A = lab[vars[idx[, 1]]], B = lab[vars[idx[, 2]]],
             cor_mean = A$mean[idx], support_pos = A$support[idx])
}))
tab$supported <- tab$support_pos >= 0.95 | tab$support_pos <= 0.05
tab$n <- nrow(dd); tab$set <- set
write.csv(tab, paste0("resultats/T9_hmsc_associations_", set, ".csv"), row.names = FALSE)
write.csv(data.frame(set = set, n = nrow(dd), psrf_beta_median = median(psrf_beta), psrf_beta_max = max(psrf_beta)),
          paste0("resultats/T9_hmsc_convergence_", set, ".csv"), row.names = FALSE)
print(tab[tab$supported, ], digits = 3, row.names = FALSE)

# variance partitioning
vp <- computeVariancePartitioning(m)
vpt <- data.frame(component = rownames(vp$vals), vp$vals, check.names = FALSE)
names(vpt)[-1] <- lab[vars]
write.csv(vpt, paste0("resultats/T9_hmsc_varpart_", set, ".csv"), row.names = FALSE)
}
