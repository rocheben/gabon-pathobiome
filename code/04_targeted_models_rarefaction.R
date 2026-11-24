#######################################################################
# 04 - Targeted models of specific associations with
#      exposure covariates, and ASA rarefaction robustness
# Manuscript items: Table S5; Figure S1 (figS1_rarefaction.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(lme4)})
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab
source_asa <- new.env(); sys.source("code/asa_functions.R", envir = source_asa)
ctrl <- glmerControl(optimizer = "bobyqa")
tidy <- function(m, name) {
  cf <- summary(m)$coefficients
  data.frame(model = name, term = rownames(cf), OR = exp(cf[, 1]),
             lo = exp(cf[, 1] - 1.96 * cf[, 2]), hi = exp(cf[, 1] + 1.96 * cf[, 2]),
             p = cf[, 4], n = nobs(m), row.names = NULL)
}
out <- list()

# ---- (10a) helminth co-infection and number of viral exposures
vir <- c("HCV","paHTLV","ebola_IgG","rift","dengue","chik")   # WNV excluded (half missing)
dv <- d[complete.cases(d[, c(vir, "Loaloa","Fil_perstans","agec","sex01","village")]), ]
dv$nvir <- rowSums(dv[, vir])
dv$helm <- factor(with(dv, ifelse(Loaloa == 1 & Fil_perstans == 1, "Both",
                           ifelse(Loaloa == 1, "L. loa only", ifelse(Fil_perstans == 1, "M. perstans only", "None")))),
                  c("None","L. loa only","M. perstans only","Both"))
m <- glmer(nvir ~ helm + agec + sex01 + (1|village), family = poisson, data = dv, control = ctrl)
out$nvir <- transform(tidy(m, "N viral exposures ~ helminth status (Poisson, rate ratio)"))
m2 <- glmer(nvir ~ Loaloa * Fil_perstans + agec + sex01 + (1|village), family = poisson, data = dv, control = ctrl)
out$nvir_int <- tidy(m2, "N viral exposures ~ Loa x Mp interaction (Poisson)")
desc_nvir <- aggregate(nvir ~ helm, dv, function(z) c(n = length(z), mean = mean(z)))
write.csv(do.call(data.frame, desc_nvir), "resultats/T4_nvir_by_helminth_descr.csv", row.names = FALSE)
# each virus separately
for (v in vir) {
  mv <- try(glmer(as.formula(paste(v, "~ Loaloa * Fil_perstans + agec + sex01 + (1|village)")),
                  family = binomial, data = dv, control = ctrl), silent = TRUE)
  if (!inherits(mv, "try-error")) out[[paste0("vir_", v)]] <- tidy(mv, paste(lab[v], "~ Loa x Mp"))
}

# ---- (9b / 10b-c) malaria models with behavioural exposure covariates
dx <- d; dx$hunter[is.na(dx$hunter)] <- NA
for (y in c("ebola_IgG","chik","paHTLV")) {
  f0 <- paste(y, "~ paMalaria + agec + sex01 + (1|village)")
  f1 <- paste(y, "~ paMalaria + agec + sex01 + hunter + eat_monkey + eat_bat + pygmee + (1|village)")
  dd <- dx[complete.cases(dx[, c(y, "paMalaria","agec","sex01","hunter","eat_monkey","eat_bat","pygmee","village")]), ]
  out[[paste0(y, "_base")]] <- tidy(glmer(as.formula(f0), family = binomial, data = dd, control = ctrl), paste(lab[y], "~ malaria (base)"))
  out[[paste0(y, "_behav")]] <- tidy(glmer(as.formula(f1), family = binomial, data = dd, control = ctrl), paste(lab[y], "~ malaria + exposure behaviour"))
}
# L. loa - HTLV with behaviour
dd <- dx[complete.cases(dx[, c("paHTLV","Loaloa","Fil_perstans","paMalaria","agec","sex01","hunter","pygmee","village")]), ]
out$htlv_loa <- tidy(glmer(paHTLV ~ Loaloa + Fil_perstans + paMalaria + agec + sex01 + hunter + pygmee + (1|village),
                           family = binomial, data = dd, control = ctrl), "HTLV ~ L. loa + M. perstans + malaria + behaviour")
# malaria parasite species vs Ebola (which species drives it)
dd <- d[complete.cases(d[, c("ebola_IgG","Pfalci","Pmalariae","Povale","agec","sex01","village")]), ]
out$ebo_sp <- tidy(glmer(ebola_IgG ~ Pfalci + Pmalariae + Povale + agec + sex01 + (1|village),
                         family = binomial, data = dd, control = ctrl), "Ebola IgG ~ Plasmodium species")
# sickle-cell trait as modifier of malaria-Ebola link
dd <- d[complete.cases(d[, c("ebola_IgG","paMalaria","AS","agec","sex01","village")]), ]
out$ebo_AS <- tidy(glmer(ebola_IgG ~ paMalaria * AS + agec + sex01 + (1|village),
                         family = binomial, data = dd, control = ctrl), "Ebola IgG ~ malaria x sickle-cell trait")

T4 <- do.call(rbind, out); T4 <- T4[T4$term != "(Intercept)", ]
write.csv(T4, "resultats/T4_targeted_models.csv", row.names = FALSE)
print(T4, digits = 3, row.names = FALSE)

# ---- ASA rarefaction for key sets
sets <- list(S2 = c("Loaloa","Fil_perstans","paHTLV"), S3 = c("Pfalci","Pmalariae","Povale"),
             S4 = c("paMalaria","dengue","westnile"), S5 = c("paMalaria","ebola_IgG","chik"))
fr <- seq(0.1, 1, 0.1); reps <- 20
rr <- list(); pb <- startpb(0, length(sets) * 2 * length(fr) * reps)
for (s in names(sets)) for (nl in c("global","village")) for (f in fr) for (r in 1:reps) {
  dd <- d[complete.cases(d[, c(sets[[s]], "village","agecl")]), ]
  set.seed(1000 * r + round(100 * f))
  sub <- dd[sample.int(nrow(dd), round(f * nrow(dd))), ]
  a <- source_asa$asa(sets[[s]], nl, data = sub, ns = 500, seed = r)
  rr[[length(rr) + 1]] <- data.frame(set = s, null = nl, frac = f, rep = r, combo = a$combo,
                                     p = a$p, sig = a$signif_envelope)
  setpb(pb, length(rr))
}
closepb(pb)
rr <- do.call(rbind, rr)
rs <- aggregate(cbind(p05 = p < 0.05, env = sig) ~ set + null + frac + combo, rr, mean)
write.csv(rs, "resultats/T5_ASA_rarefaction.csv", row.names = FALSE)
# minimum sampling fraction from which >=80% of replicates are significant at every larger fraction
minfrac <- do.call(rbind, lapply(split(rs, list(rs$set, rs$null, rs$combo), drop = TRUE), function(z) {
  z <- z[order(z$frac), ]; ok <- z$p05 >= 0.8
  mf <- if (!ok[length(ok)]) NA else { i <- max(which(!ok), 0) + 1; z$frac[i] }
  data.frame(set = z$set[1], null = z$null[1], combo = z$combo[1], p05_full = z$p05[nrow(z)], min_frac_80 = mf)
}))
write.csv(minfrac, "resultats/T5_ASA_rarefaction_minfrac.csv", row.names = FALSE)
