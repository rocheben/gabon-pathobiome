#######################################################################
# 12 - Ecological scales of co-occurrence (tables T14, T15)
#  (a) Hmsc variance partitioning (9 pathogen groups): share of explained
#      variation attributed to ecosystem, host (age, sex), village-level
#      and individual-level latent factors
#  (b) For all pairs of the 10 pathogen groups: village-scale association
#      (Spearman correlation of village prevalences, partial on ecosystem)
#      versus individual-scale association (adjusted log OR, age, sex and
#      village). Concordance between scales.
# Manuscript items: Table S8; Figure 2 (fig2_scales.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
P <- readRDS("resultats/prepped.rds"); lab <- P$lab
vp <- read.csv("resultats/T9_hmsc_varpart_9.csv", check.names = FALSE)
comp <- c(age = "Host (age, sex)", `I(age^2)` = "Host (age, sex)", sex = "Host (age, sex)",
          eco = "Ecosystem", `Random: village` = "Village (latent)", `Random: indiv` = "Individual (latent)")
long <- do.call(rbind, lapply(names(vp)[-1], function(s)
  data.frame(pathogen = s, component = comp[vp$component], share = vp[[s]])))
long <- aggregate(share ~ pathogen + component, long, sum); long$pathogen <- sub("Ebola (IgG)", "EBOV IgG", long$pathogen, fixed = TRUE)
write.csv(long, "resultats/T14_hmsc_varpart_summary.csv", row.names = FALSE)
print(aggregate(share ~ component, long, mean))
# (b) village vs individual scale
vprev <- read.csv("resultats/T7_village_prevalence.csv")
pw <- read.csv("resultats/T2_pairwise.csv")
p10 <- P$patho10
pairs <- t(combn(p10, 2))
sc <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  z <- vprev[complete.cases(vprev[, c(a, b)]), ]
  if (nrow(z) < 20) return(NULL)
  r1 <- resid(lm(rank(z[[a]]) ~ z$ecosys2)); r2 <- resid(lm(rank(z[[b]]) ~ z$ecosys2))
  ct <- cor.test(r1, r2)
  row <- pw[(pw$A == a & pw$B == b) | (pw$A == b & pw$B == a), ]
  data.frame(A = lab[a], B = lab[b], n_villages = nrow(z), rho_village = unname(ct$estimate),
             p_village = ct$p.value, logOR_individual = log(max(row$OR_adj, 0.05)), q_individual = row$q_adj)
}))
sc$q_village <- p.adjust(sc$p_village, "BH")
sc$class <- with(sc, ifelse(q_village < 0.05 & q_individual < 0.05, "Both scales",
                     ifelse(q_village < 0.05, "Village scale only",
                     ifelse(q_individual < 0.05, "Individual scale only", "Neither"))))
write.csv(sc, "resultats/T15_scale_concordance.csv", row.names = FALSE)
ct <- cor.test(sc$rho_village, sc$logOR_individual, method = "spearman")
write.csv(data.frame(n_pairs = nrow(sc), rho = unname(ct$estimate), p = ct$p.value,
                     table(sc$class)), "resultats/T15_scale_concordance_summary.csv", row.names = FALSE)
print(table(sc$class)); print(ct)
