#######################################################################
# Figure S1 - Rarefaction of association screening (resultats/T5_ASA_rarefaction*.csv)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
P <- readRDS("resultats/prepped.rds"); lab <- P$lab
sets <- list(S2 = c("Loaloa","Fil_perstans","paHTLV"), S3 = c("Pfalci","Pmalariae","Povale"),
             S4 = c("paMalaria","dengue","westnile"), S5 = c("paMalaria","ebola_IgG","chik"))
rs <- read.csv("resultats/T5_ASA_rarefaction.csv"); minfrac <- read.csv("resultats/T5_ASA_rarefaction_minfrac.csv")
rs$setlab <- sapply(sets[rs$set], function(v) paste(lab[v], collapse = " / "))
rs$combo <- sub("^c", "", rs$combo)
g <- ggplot(rs[paste(rs$set, rs$combo) %in% paste(minfrac$set, sub("^c", "", minfrac$combo))[!is.na(minfrac$min_frac_80)], ],
            aes(frac, p05, colour = combo, linetype = null)) + geom_line() +
  facet_wrap(~setlab) + geom_hline(yintercept = 0.8, linetype = 3) +
  labs(x = "Fraction of hosts subsampled", y = "Share of subsamples with P < 0.05", colour = "Combination", linetype = "Null model") +
  theme_bw(base_size = 8) + theme(strip.text = element_text(face = "italic"))
save_fig("FigS1_rarefaction", function() print(g), width = 6.4, height = 5)
