#######################################################################
# 05 - Bayesian network analysis
#  - bootstrap arc strength (500 resamples, hill-climbing, BDe iss=1)
#  - model-averaged network (Scutari & Nagarajan 2013 threshold)
#  - same with environmental/host nodes (province, ecosystem, age class,
#    sex): arcs between pathogens that persist are not explained by those
#  - sign of each arc from the conditional odds ratio
# Manuscript items: Appendix S1 Section S6; Figure S4 (figS4_bayesian_networks.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(bnlearn)})
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab
set.seed(42); R <- 500

run_bn <- function(vars, covs = NULL, name) {
  dd <- d[complete.cases(d[, c(vars, covs)]), c(vars, covs)]
  dd[] <- lapply(dd, factor)
  bl <- NULL
  if (!is.null(covs)) bl <- expand.grid(from = vars, to = covs, stringsAsFactors = FALSE)
  bl <- rbind(bl, if (length(covs) > 1) {x <- expand.grid(from = covs, to = covs, stringsAsFactors = FALSE); x[x$from != x$to, ]})
  bs <- boot.strength(dd, R = R, algorithm = "hc",
                      algorithm.args = list(score = "bde", iss = 1, blacklist = bl))
  thr <- attr(bs, "threshold")
  und <- bs[bs$from < bs$to, ]
  # direction-free strength = strength; keep pathogen-pathogen arcs
  und <- und[und$from %in% vars & und$to %in% vars, ]
  und$sign <- mapply(function(a, b) {
    t <- table(dd[[a]], dd[[b]]); (t[2,2] * t[1,1]) / (t[1,2] * t[2,1])}, und$from, und$to)
  und$set <- name; und$n <- nrow(dd); und$threshold <- thr
  und$in_avg <- und$strength >= thr
  list(tab = und, avg = averaged.network(bs), n = nrow(dd), dd = dd, bl = bl)
}

p10 <- P$patho10; p12 <- P$patho12; p9 <- setdiff(p10, "westnile")
covs <- c("Province", "ecosys2", "agecl", "sex")
specs <- list(
  bn10       = list(p10, NULL, "10 taxa"),
  bn10_env   = list(p10, covs, "10 taxa + environment/host"),
  bn12       = list(p12, NULL, "12 taxa"),
  bn12_env   = list(p12, covs, "12 taxa + environment/host"),
  bn9_env    = list(p9,  covs, "9 taxa (no WNV) + environment/host")
)
runs <- pblapply(specs, function(s) run_bn(s[[1]], s[[2]], s[[3]]))
tab <- do.call(rbind, lapply(runs, `[[`, "tab"))
tab$labFrom <- lab[tab$from]; tab$labTo <- lab[tab$to]
tab <- tab[tab$strength >= 0.3, ]
tab <- tab[order(tab$set, -tab$strength), ]
write.csv(tab, "resultats/T6_BN_arc_strength.csv", row.names = FALSE)
print(tab[tab$in_avg, c("set","n","labFrom","labTo","strength","sign","threshold")], digits = 3, row.names = FALSE)

