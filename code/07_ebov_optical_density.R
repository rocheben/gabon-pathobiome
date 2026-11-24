#######################################################################
# 07 - Ebola IgG optical density (ebola_DO) by malaria status:
#      does malaria shift the whole OD distribution (non-specific
#      reactivity) or only the proportion of clear positives?
# Manuscript items: Table S9; Figure S7 (figS7_EBOV_optical_density.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(lme4))
P <- readRDS("resultats/prepped.rds"); d <- P$d
d$ebola_DO <- suppressWarnings(as.numeric(d$ebola_DO))
z <- d[complete.cases(d[, c("ebola_DO","ebola_IgG","paMalaria","agec","sex01","village")]), ]
s <- aggregate(ebola_DO ~ paMalaria + ebola_IgG, z, function(x) c(n = length(x), median = median(x), q1 = quantile(x, .25), q3 = quantile(x, .75)))
s <- do.call(data.frame, s); print(s, digits = 3)
w_neg <- wilcox.test(ebola_DO ~ paMalaria, z[z$ebola_IgG == 0, ])
w_pos <- wilcox.test(ebola_DO ~ paMalaria, z[z$ebola_IgG == 1, ])
m <- lmer(log(ebola_DO) ~ paMalaria + agec + sex01 + (1|village), data = z[z$ebola_IgG == 0 & z$ebola_DO > 0, ])
cf <- summary(m)$coefficients
thr <- range(z$ebola_DO[z$ebola_IgG == 1]); thr0 <- range(z$ebola_DO[z$ebola_IgG == 0])
res <- data.frame(test = c("Wilcoxon OD malaria+ vs - among IgG negatives", "Wilcoxon OD malaria+ vs - among IgG positives",
                           "LMM log(OD) ~ malaria among IgG negatives (ratio)"),
                  estimate = c(NA, NA, exp(cf["paMalaria", 1])), p = c(w_neg$p.value, w_pos$p.value,
                  2 * pnorm(-abs(cf["paMalaria", 3]))))
print(res, digits = 3); cat("OD range positives:", thr, " negatives:", thr0, "\n")
write.csv(s, "resultats/T8_ebolaOD_by_malaria.csv", row.names = FALSE)
write.csv(res, "resultats/T8_ebolaOD_tests.csv", row.names = FALSE)

# sensitivity: malaria-Ebola adjusted OR under stricter OD cut-offs
sens <- do.call(rbind, lapply(c(0.2, 0.3, 0.5, 0.8, 1.0), function(ct) {
  z$eb <- as.integer(z$ebola_DO >= ct)
  m <- glmer(eb ~ paMalaria + agec + sex01 + (1|village), family = binomial, data = z,
             control = glmerControl(optimizer = "bobyqa"))
  cf <- summary(m)$coefficients["paMalaria", ]
  data.frame(cutoff = ct, n_pos = sum(z$eb), OR = exp(cf[1]), lo = exp(cf[1] - 1.96 * cf[2]),
             hi = exp(cf[1] + 1.96 * cf[2]), p = cf[4])
}))
print(sens, digits = 3, row.names = FALSE)
write.csv(sens, "resultats/T8_ebola_cutoff_sensitivity.csv", row.names = FALSE)
