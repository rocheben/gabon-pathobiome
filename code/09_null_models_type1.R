#######################################################################
# 09 - Association screening: FDR across combinations and type I error of the null models
#  A. All 75 pairs refitted with penalised GLMMs (blme, weakly informative
#     priors; handles separation) and a natural-spline age effect
#  B. Key associations: extended exposure covariates, conditioning on the
#     rest of the pathogen community, exclusion of lagoon villages,
#     leave-one-province-out, E-values for unmeasured confounding
#  C. False-discovery control across combinations in association screening
#  D. Simulation: false-positive rate of the three ASA null models when
#     pathogens are independent within individuals but prevalence varies
#     among villages and age classes (the confounding scenario)
# Manuscript items: Table S3 (FDR); Table S4; Figure S2a (figS2_null_benchmark.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(blme); library(splines)})
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab
source("code/asa_functions.R", local = TRUE)
d$province <- factor(d$Province)
ctrl <- glmerControl(optimizer = "bobyqa")

bfit <- function(y, x, dd, extra = "") {
  f <- as.formula(paste(y, "~", x, "+ ns(age, 3) + sex01", extra, "+ (1|village)"))
  m <- try(bglmer(f, data = dd, family = binomial, control = ctrl,
                  fixef.prior = normal(sd = c(10, 2.5)), cov.prior = gamma), silent = TRUE)
  if (inherits(m, "try-error")) return(c(est = NA, se = NA, p = NA))
  cf <- summary(m)$coefficients
  c(est = cf[x, 1], se = cf[x, 2], p = cf[x, 4])
}
sym <- function(a, b, dd, extra = "") {
  f1 <- bfit(a, b, dd, extra); f2 <- bfit(b, a, dd, extra)
  w1 <- 1 / f1["se"]^2; w2 <- 1 / f2["se"]^2
  est <- (f1["est"] * w1 + f2["est"] * w2) / (w1 + w2); se <- sqrt(1 / (w1 + w2))
  c(OR = unname(exp(est)), lo = unname(exp(est - 1.96 * se)), hi = unname(exp(est + 1.96 * se)),
    p = unname(max(f1["p"], f2["p"])))
}

# ---- C. BH across combinations within set x null
asa_tab <- read.csv("resultats/T3_ASA_all.csv")
asa_tab$q <- ave(asa_tab$p, asa_tab$set_id, asa_tab$null, FUN = function(z) p.adjust(z, "BH"))
asa_tab$robust_all3 <- ave(asa_tab$signif_envelope & asa_tab$q < 0.05, asa_tab$set_id, asa_tab$combo, FUN = all)
write.csv(asa_tab, "resultats/T3_ASA_all_withFDR.csv", row.names = FALSE)

# ---- D. simulation of ASA false-positive rates under heterogeneous exposure
sim_sets <- list(S1 = c("Loaloa","Fil_perstans","paMalaria"), S5 = c("paMalaria","ebola_IgG","chik"))
nrep <- 50
D <- do.call(rbind, lapply(names(sim_sets), function(s) {
  v <- sim_sets[[s]]
  dd <- d[complete.cases(d[, c(v, "village", "agecl")]), ]
  st <- interaction(dd$village, dd$agecl, drop = TRUE)
  pr <- sapply(v, function(z) ave(dd[[z]], st))          # stratum prevalence
  message("Simulation, set ", s)
  do.call(rbind, pblapply(1:nrep, function(r) {
    set.seed(500 + r); sim <- dd
    for (j in seq_along(v)) sim[[v[j]]] <- rbinom(nrow(dd), 1, pr[, j])
    do.call(rbind, lapply(c("global","village","vill_age"), function(nl) {
      a <- asa(v, nl, data = sim, ns = 500, seed = r)
      data.frame(set = s, rep = r, null = nl, any_sig = any(a$signif_envelope),
                 any_p05 = any(p.adjust(a$p, "BH") < 0.05))
    }))
  }))
}))
Dsum <- aggregate(cbind(any_sig, any_p05) ~ set + null, D, mean)
names(Dsum)[3:4] <- c("false_positive_rate_envelope", "false_positive_rate_BH")
write.csv(Dsum, "resultats/T12_ASA_null_simulation.csv", row.names = FALSE)
print(Dsum)
