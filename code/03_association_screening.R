#######################################################################
# 03 - Association screening analysis (ASA; Vaumourin et al. 2014,
#      Gasqui FctTestScreenENV) extended with exposure-controlled nulls
#
#  Null "global"   : Gasqui's original null. Combination counts drawn from a
#                    multinomial with probabilities = product of global
#                    prevalences (independence, homogeneous population).
#  Null "village"  : each pathogen permuted independently WITHIN village.
#                    Keeps every village's prevalence, removes shared-exposure
#                    signal due to geography/habitat.
#  Null "vill_age" : permutation within village x age class (<=30, 31-45,
#                    46-60, >60). Also removes cumulative-exposure (age) signal
#                    of serological markers.
#  95% global confidence envelope (boot::envelope) as in the original code.
# Manuscript items: Table S3; Figure 4 (fig4_association_screening.R)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages({library(boot)})
source("code/setup.R")  # installs missing packages, loads pbapply (progress bars)
P <- readRDS("resultats/prepped.rds"); d <- P$d; lab <- P$lab


source("code/asa_functions.R", local = TRUE)

sets <- list(
  S1_helminths_malaria   = c("Loaloa","Fil_perstans","paMalaria"),
  S2_helminths_HTLV      = c("Loaloa","Fil_perstans","paHTLV"),
  S3_plasmodium          = c("Pfalci","Pmalariae","Povale"),
  S4_malaria_DENV_WNV    = c("paMalaria","dengue","westnile"),
  S5_malaria_ebola_chik  = c("paMalaria","ebola_IgG","chik"),
  S6_helminths_viruses   = c("Loaloa","Fil_perstans","HCV","paHTLV","ebola_IgG"),
  S7_arboviruses         = c("rift","dengue","westnile","chik"),
  S8_all10               = P$patho10
)
res <- list(); pb <- startpb(0, 3 * length(sets))
for (s in names(sets)) for (nl in c("global","village","vill_age")) {
  r <- asa(sets[[s]], nl); r$set_id <- s; res[[paste(s, nl)]] <- r; setpb(pb, length(res))
}
closepb(pb)
res <- do.call(rbind, res)
write.csv(res, "resultats/T3_ASA_all.csv", row.names = FALSE)

# summary table: combinations significant under at least one null
w <- reshape(res[, c("set_id","combo","N","obs","null","expected","direction","p")],
             idvar = c("set_id","combo"), timevar = "null", direction = "wide", drop = "N")
keep <- with(w, direction.global != "random" | direction.village != "random" | direction.vill_age != "random")
write.csv(w[keep, ], "resultats/T3_ASA_signif_summary.csv", row.names = FALSE)
print(w[keep, ], digits = 3, row.names = FALSE)

