#######################################################################
# 10 - Benchmark of association-screening null models: statistical power
#  Data are simulated with the observed village-by-age prevalences of
#  Plasmodium, EBOV IgG and CHIKV (set S5), and a within-host dependence
#  of known strength (odds ratio) is injected between the first two.
#  Power = share of simulated datasets in which the A+B co-infection
#  combinations (110 and 111) are flagged in excess outside the envelope.
#  Type I error: see 09_null_models_type1.R (T12).
# Manuscript items: Table S4; Figure S2b (figS2_null_benchmark.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)
source("code/asa_functions.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); d <- P$d
v <- c("paMalaria", "ebola_IgG", "chik")
dd <- d[complete.cases(d[, c(v, "village", "agecl")]), ]
st <- interaction(dd$village, dd$agecl, drop = TRUE)
pr <- sapply(v, function(z) ave(dd[[z]], st))
clamp <- function(p) pmin(pmax(p, 1e-4), 1 - 1e-4)
nrep <- 30
res <- do.call(rbind, lapply(c(1.5, 2, 3), function(or) {
  message("Power simulation, OR = ", or)
  do.call(rbind, pblapply(1:nrep, function(r) {
    set.seed(9000 + 100 * or + r); sim <- dd
    A <- rbinom(nrow(dd), 1, pr[, 1])
    B <- rbinom(nrow(dd), 1, plogis(qlogis(clamp(pr[, 2])) + log(or) * (A - pr[, 1])))
    C <- rbinom(nrow(dd), 1, pr[, 3])
    sim$paMalaria <- A; sim$ebola_IgG <- B; sim$chik <- C
    do.call(rbind, lapply(c("global", "village", "vill_age"), function(nl) {
      a <- asa(v, nl, data = sim, ns = 500, seed = r)
      hit <- a$combo %in% c("c110", "c111") & a$direction == "frequent"
      data.frame(OR = or, rep = r, null = nl, detected = any(hit))
    }))
  }))
}))
pw <- aggregate(detected ~ OR + null, res, mean)
names(pw)[3] <- "power"
write.csv(pw, "resultats/T13_ASA_null_power.csv", row.names = FALSE)
print(pw)
