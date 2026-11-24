#######################################################################
# 08 - Robustness analyses
#  A. All 75 pairs refitted with penalised GLMMs (blme, weakly informative
#     priors; handles separation) and a natural-spline age effect
#  B. Key associations: extended exposure covariates, conditioning on the
#     rest of the pathogen community, exclusion of lagoon villages,
#     leave-one-province-out, E-values for unmeasured confounding
#  C. False-discovery control across combinations in association screening
#  D. Simulation: false-positive rate of the three ASA null models when
#     pathogens are independent within individuals but prevalence varies
#     among villages and age classes (the confounding scenario)
# Manuscript items: Table 1; Table S2 (penalised models); Figure S3 (figS3_robustness.R)
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

# ---- A. all pairs, penalised GLMM + spline age
prim <- read.csv("resultats/T2_pairwise.csv")
A <- do.call(rbind, pblapply(seq_len(nrow(prim)), function(i) {
  a <- prim$A[i]; b <- prim$B[i]
  dd <- d[complete.cases(d[, c(a, b, "age", "sex01", "village")]), ]
  r <- sym(a, b, dd); data.frame(A = a, B = b, t(r))
}))
A$q <- p.adjust(A$p, "BH")
A <- merge(prim[, c("A","B","labA","labB","n","OR_adj","q_adj","status")], A, by = c("A","B"))
A$agree <- (A$q < 0.05) == (A$q_adj < 0.05)
write.csv(A, "resultats/T10_pairs_penalised_spline.csv", row.names = FALSE)
cat("Agreement primary vs penalised/spline:", sum(A$agree), "/", nrow(A), "\n")
print(A[!A$agree | A$q < 0.05, c("labA","labB","OR_adj","q_adj","OR","lo","hi","q")], digits = 3)

# ---- B. key associations
key <- list(c("Pmalariae","Pfalci"), c("Povale","Pmalariae"), c("westnile","dengue"),
            c("ebola_IgG","paMalaria"), c("ebola_IgG","Pfalci"), c("ebola_IgG","Pmalariae"),
            c("paHTLV","paMalaria"), c("paHTLV","Pfalci"), c("paHTLV","Loaloa"),
            c("Fil_perstans","Loaloa"))
others9 <- c("Loaloa","Fil_perstans","paMalaria","HCV","paHTLV","ebola_IgG","rift","dengue","chik")
behav <- "+ hunter + pygmee + eat_monkey + eat_bat"
plas <- c("paMalaria","Pfalci","Pmalariae","Povale")
evalue <- function(or, common) { rr <- if (common) sqrt(or) else or
  if (rr < 1) rr <- 1 / rr; rr + sqrt(rr * (rr - 1)) }
B <- do.call(rbind, pblapply(key, function(pr) {
  y <- pr[1]; x <- pr[2]
  base <- d[complete.cases(d[, c(y, x, "age", "sex01", "village")]), ]
  r_sym <- sym(y, x, base)
  # extended behaviour
  db <- base[complete.cases(base[, c("hunter","pygmee","eat_monkey","eat_bat")]), ]
  r_beh <- bfit(y, x, db, behav)
  # conditional on the rest of the community (other pathogens as covariates)
  oth <- setdiff(others9, c(y, x, if (x %in% plas | y %in% plas) plas))
  if (y == "westnile" | x == "westnile") oth <- setdiff(oth, "dengue")
  dc <- base[complete.cases(base[, oth]), ]
  r_con <- bfit(y, x, dc, paste(c("", oth), collapse = " + "))
  # without lagoon
  dl <- base[base$ecosys2 != "LAGUNE", ]; r_lag <- bfit(y, x, dl)
  # leave one province out
  lopo <- sapply(levels(droplevels(base$province)), function(pv) {
    r <- bfit(y, x, base[base$province != pv, ]); c(exp(r["est"]), r["p"]) })
  common <- mean(base[[y]]) > 0.15 | mean(base[[x]]) > 0.15
  lo_ci <- r_sym["lo"]; e_pt <- evalue(r_sym["OR"], common)
  e_ci <- if ((r_sym["lo"] > 1) | (r_sym["hi"] < 1)) evalue(ifelse(r_sym["OR"] > 1, r_sym["lo"], r_sym["hi"]), common) else 1
  data.frame(outcome = lab[y], exposure = lab[x], n = nrow(base),
             OR = r_sym["OR"], lo = r_sym["lo"], hi = r_sym["hi"], p = r_sym["p"],
             OR_behaviour = exp(r_beh["est"]), p_behaviour = r_beh["p"], n_behaviour = nrow(db),
             OR_conditional = exp(r_con["est"]), p_conditional = r_con["p"], n_conditional = nrow(dc),
             OR_no_lagoon = exp(r_lag["est"]), p_no_lagoon = r_lag["p"],
             LOPO_OR_min = min(lopo[1, ], na.rm = TRUE), LOPO_OR_max = max(lopo[1, ], na.rm = TRUE),
             LOPO_n_sig = sum(lopo[2, ] < 0.05, na.rm = TRUE), LOPO_n = ncol(lopo),
             Evalue_point = e_pt, Evalue_CI = e_ci, row.names = NULL)
}))
write.csv(B, "resultats/T11_key_associations_robustness.csv", row.names = FALSE)
print(B, digits = 3)

