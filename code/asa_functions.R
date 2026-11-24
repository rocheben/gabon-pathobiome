# Association screening functions (used by 03, 04, 09, 10)
suppressMessages(library(boot))
NS <- 5000
combo_code <- function(X) as.vector(X %*% 2^((ncol(X) - 1):0))
combo_labels <- function(vars) {
  k <- length(vars); sapply(0:(2^k - 1), function(z) paste(rev(as.integer(intToBits(z))[1:k]), collapse = ""))
}

asa <- function(vars, null = c("global", "village", "vill_age"), data = d, ns = NS, seed = 1) {
  null <- match.arg(null)
  dd <- data[complete.cases(data[, c(vars, "village", "agecl")]), ]
  strata <- switch(null, global = rep(1, nrow(dd)), village = as.integer(dd$village),
                   vill_age = as.integer(interaction(dd$village, dd$agecl, drop = TRUE)))
  o <- order(strata); dd <- dd[o, ]; strata <- strata[o]
  X <- as.matrix(dd[, vars]); N <- nrow(X); k <- length(vars); nc <- 2^k
  obs <- tabulate(combo_code(X) + 1, nbins = nc)
  set.seed(seed)
  if (null == "global") {
    pv <- colMeans(X)
    B <- as.matrix(expand.grid(rep(list(1:0), k)))[, k:1, drop = FALSE]  # rows = combos in code order
    B <- B[order(combo_code(B)), , drop = FALSE]
    pcomb <- apply(B, 1, function(b) prod(ifelse(b == 1, pv, 1 - pv)))
    sims <- t(rmultinom(ns, N, pcomb))
  } else {
    sims <- t(vapply(seq_len(ns), function(s) {
      Xp <- X
      for (j in seq_len(k)) Xp[, j] <- X[order(strata, runif(N)), j]
      tabulate(combo_code(Xp) + 1, nbins = nc)
    }, numeric(nc)))
  }
  env <- envelope(mat = sims, level = 0.95)
  hi <- env$overall[1, ]; lo <- env$overall[2, ]
  exp_mean <- colMeans(sims)
  p_hi <- (colSums(sims >= rep(obs, each = ns)) + 1) / (ns + 1)
  p_lo <- (colSums(sims <= rep(obs, each = ns)) + 1) / (ns + 1)
  pval <- pmin(1, 2 * pmin(p_hi, p_lo))
  sig <- obs > hi | obs < lo
  data.frame(set = paste(vars, collapse = "+"), null = null, N = N,
             combo = paste0("c", combo_labels(vars)), obs = obs, expected = exp_mean,
             env_lo = lo, env_hi = hi, ratio = obs / exp_mean, p = pval,
             signif_envelope = sig,
             direction = ifelse(!sig, "random", ifelse(obs > hi, "frequent", "rare")))
}

