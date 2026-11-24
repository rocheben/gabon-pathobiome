#######################################################################
# Figure S6 - Village prevalence of Plasmodium and EBOV IgG (resultats/T7_village_prevalence.csv)
# Run from the project root (folder containing code/, data/, resultats/, figures/), e.g. source("code/run_all.R").
#######################################################################
suppressMessages(library(ggplot2)); source("code/save_figure.R", local = TRUE)
vp <- read.csv("resultats/T7_village_prevalence.csv")
z <- vp[complete.cases(vp[, c("ebola_IgG","paMalaria")]), ]
g <- ggplot(z, aes(paMalaria, ebola_IgG, colour = ecosys2, size = n)) + geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(FORET = "#1b7837", SAVANE = "#d9a520", LAGUNE = "#2166ac"),
                      breaks = c("FORET","SAVANE","LAGUNE"), labels = c("Forest","Savanna","Lagoon"), name = NULL) +
  scale_size_area(max_size = 4, name = "Hosts sampled") +
  labs(x = "Plasmodium prevalence (village)", y = "EBOV IgG seroprevalence (village)") +
  theme_classic(base_size = 9)
save_fig("FigS6_village_EBOV_Plasmodium", function() print(g), width = 5.5, height = 4)
