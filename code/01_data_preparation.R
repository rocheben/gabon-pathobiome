#######################################################################
# 01 - Data preparation and descriptive analysis of the pathogen community
# Gabon pathobiome (Abbate et al.)
# Manuscript items: Table S1; Figure 1 (tables read by fig1_community_structure.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
dir.create("resultats", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)

raw <- read.csv("data/gabondata_master_utf8.csv", stringsAsFactors = FALSE)

# Pathogen variables (12 taxa; malaria also as pooled Plasmodium)
patho12 <- c("Loaloa","Fil_perstans","Pfalci","Pmalariae","Povale",
             "HCV","paHTLV","ebola_IgG","rift","dengue","westnile","chik")
patho10 <- c("Loaloa","Fil_perstans","paMalaria","HCV","paHTLV",
             "ebola_IgG","rift","dengue","westnile","chik")
lab <- c(Loaloa="L. loa", Fil_perstans="M. perstans", paMalaria="Plasmodium spp.",
         Pfalci="P. falciparum", Pmalariae="P. malariae", Povale="P. ovale",
         HCV="HCV", paHTLV="HTLV", ebola_IgG="EBOV IgG", rift="RVFV",
         dengue="DENV", westnile="WNV", chik="CHIKV")
guild <- c(Loaloa="Helminth", Fil_perstans="Helminth", paMalaria="Protozoa",
           Pfalci="Protozoa", Pmalariae="Protozoa", Povale="Protozoa",
           HCV="Virus (blood-borne)", paHTLV="Virus (blood-borne)", ebola_IgG="Virus (zoonotic)",
           rift="Virus (arbovirus)", dengue="Virus (arbovirus)", westnile="Virus (arbovirus)",
           chik="Virus (arbovirus)")
d <- raw
d$sex01   <- as.integer(tolower(d$sex) == "masculin")
d$agec    <- as.numeric(scale(d$age))
d$agecl   <- cut(d$age, c(0, 30, 45, 60, 120), labels = c("<=30","31-45","46-60",">60"))
d$ecosys2 <- factor(d$ecosys2)
d$village <- factor(d$village)
d$dept    <- factor(d$dept)
saveRDS(list(d = d, patho12 = patho12, patho10 = patho10, lab = lab, guild = guild), "resultats/prepped.rds")

# ---- Table S1: prevalence per pathogen
tab1 <- do.call(rbind, lapply(c(patho10, "Pfalci","Pmalariae","Povale"), function(p) {
  x <- d[[p]]; n <- sum(!is.na(x)); k <- sum(x == 1, na.rm = TRUE)
  ci <- binom.test(k, n)$conf.int
  data.frame(var = p, label = lab[p], guild = guild[p], n = n, pos = k,
             prev = k / n, lo = ci[1], hi = ci[2],
             n_villages = length(unique(d$village[!is.na(x)])))
}))
write.csv(tab1, "resultats/T1_prevalence.csv", row.names = FALSE)

# ---- consistency check with the pre-built complete-case files
cc10 <- complete.cases(d[, patho10]); cc12 <- complete.cases(d[, patho12])
chk <- data.frame(set = c("Gabon10All", "Gabon12All", "Gabon10Depts"),
                  n_file = c(nrow(read.csv("data/Gabon10All.csv")), nrow(read.csv("data/Gabon12All.csv")),
                             nrow(read.csv("data/Gabon10Depts.csv"))),
                  n_master_complete = c(sum(cc10), sum(cc12), NA))
write.csv(chk, "resultats/check_completecases.csv", row.names = FALSE)

# prevalence by ecosystem (for habitat section), complete-case 10
prev_eco <- do.call(rbind, lapply(patho10, function(p) {
  tt <- table(d$ecosys2, d[[p]])
  data.frame(var = p, label = lab[p], ecosys = rownames(tt), n = rowSums(tt),
             prev = tt[, "1"] / rowSums(tt), p_chisq = suppressWarnings(chisq.test(tt)$p.value))
}))
write.csv(prev_eco, "resultats/T_prevalence_by_ecosystem.csv", row.names = FALSE)

# ---- Number of infections per individual vs independence (Poisson-binomial)
X <- as.matrix(d[cc10, patho10])
p <- colMeans(X)
obs <- table(factor(rowSums(X), levels = 0:length(p)))
# exact Poisson-binomial by convolution
pb <- 1; for (pi in p) pb <- c(pb * (1 - pi), 0) + c(0, pb * pi)
expd <- pb * nrow(X)
# null that keeps village-level prevalences (permutation within village)
set.seed(1)
vill <- droplevels(d$village[cc10])
perm_counts <- pbreplicate(2000, {
  Xp <- apply(X, 2, function(col) ave(col, vill, FUN = function(z) z[sample.int(length(z))]))
  tabulate(rowSums(Xp) + 1, nbins = length(p) + 1)
})
richness <- data.frame(k = 0:length(p), observed = as.numeric(obs), expected_indep = expd,
                       expected_village = rowMeans(perm_counts),
                       village_lo = apply(perm_counts, 1, quantile, 0.025),
                       village_hi = apply(perm_counts, 1, quantile, 0.975))
write.csv(richness, "resultats/T_richness_distribution.csv", row.names = FALSE)
vr <- c(obs_var = var(rowSums(X)), indep_var = sum(p * (1 - p)),
        village_var_mean = mean(apply(perm_counts, 2, function(cn) {
          k <- rep(0:length(p), cn); var(k)})))
write.csv(t(vr), "resultats/T_richness_variance.csv", row.names = FALSE)

cat("done\n")
