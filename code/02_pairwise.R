#######################################################################
# 02 - Pairwise associations: crude (Fisher) and adjusted (GLMM with
#      age, sex and a village random effect = shared exposure control)
# Manuscript items: Table S2; Figure 3 (fig3_pairwise.R); Figure 2b (via 12_scales.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(lme4)})
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab; guild <- P$guild
vars <- c(P$patho12, "paMalaria")
pairs <- t(combn(vars, 2))
plas <- c("paMalaria","Pfalci","Pmalariae","Povale")
pairs <- pairs[!(pairs[,1] %in% plas & pairs[,2] %in% plas & (pairs[,1]=="paMalaria" | pairs[,2]=="paMalaria")), ]

fit_one <- function(y, x, dd) {
  m <- try(glmer(as.formula(paste(y, "~", x, "+ agec + sex01 + (1|village)")), data = dd,
                 family = binomial, control = glmerControl(optimizer = "bobyqa")), silent = TRUE)
  if (inherits(m, "try-error")) return(c(NA, NA, NA, NA))
  cf <- summary(m)$coefficients
  c(cf[x, 1], cf[x, 2], cf[x, 4], as.numeric(isSingular(m)))
}

res <- do.call(rbind, pblapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  dd <- d[complete.cases(d[, c(a, b, "agec", "sex01", "village")]), ]
  tt <- table(factor(dd[[a]], 0:1), factor(dd[[b]], 0:1))
  ft <- fisher.test(tt)
  f1 <- fit_one(a, b, dd); f2 <- fit_one(b, a, dd)
  data.frame(A = a, B = b, n = nrow(dd), n11 = tt["1","1"],
             exp11 = sum(tt["1",]) * sum(tt[,"1"]) / sum(tt),
             OR_crude = unname(ft$estimate), p_crude = ft$p.value,
             logOR_adjA = f1[1], se_adjA = f1[2], p_adjA = f1[3],
             logOR_adjB = f2[1], se_adjB = f2[2], p_adjB = f2[3])
}))
# symmetric summary of the adjusted effect: inverse-variance mean of both directions
w1 <- 1 / res$se_adjA^2; w2 <- 1 / res$se_adjB^2
res$OR_adj <- exp((res$logOR_adjA * w1 + res$logOR_adjB * w2) / (w1 + w2))
res$p_adj  <- pmax(res$p_adjA, res$p_adjB)       # conservative: both directions must agree
res$q_crude <- p.adjust(res$p_crude, "BH")
res$q_adj   <- p.adjust(res$p_adj, "BH")
res$labA <- lab[res$A]; res$labB <- lab[res$B]
res$type <- ifelse(guild[res$A] == guild[res$B] |
                     (grepl("Virus", guild[res$A]) & grepl("Virus", guild[res$B])),
                   "within-guild", "between-guild")
res$status <- with(res, ifelse(q_crude < 0.05 & q_adj < 0.05, "robust",
                        ifelse(q_crude < 0.05 & q_adj >= 0.05, "explained by shared exposure",
                        ifelse(q_crude >= 0.05 & q_adj < 0.05, "revealed after adjustment", "ns"))))
res <- res[order(res$p_adj), ]
write.csv(res, "resultats/T2_pairwise.csv", row.names = FALSE)

print(res[res$status != "ns", c("labA","labB","n","n11","exp11","OR_crude","q_crude","OR_adj","q_adj","status")], digits = 3)
