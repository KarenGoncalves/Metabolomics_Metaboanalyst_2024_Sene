# MetaboAnalyst multivariate analyses
# install.packages(c("ellipse", "pls"))
library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
source("scripts/help_script/calculate.pls.vip.R")
theme_classic() %>% theme_set()

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
inFiles <- paste0("Results/", Analysis_modes, "_mSet.RData")
plot_basename <- paste0("plots/", Analysis_modes, "/")
# i = 1

for (i in 1:3) {
    rm("mSet")
    load(inFiles[i])
    mSet <- SAM.Anal(mSetObj = mSet, 
                     method = "d.stat", # d.stat for parametric 
                     paired = F,
                     varequal = T, 
                     imgName = paste0(plot_basename, 
                                      "Differential_analysis_"),
                     dpi=300)
    mSet<-PlotSAM.FDR(mSet, 
                      imgName = paste0(plot_basename, 
                                       "Differential_analysis_FDR"),
                      format = "png", 
                      dpi=300, width=NA)
    ?PlotSAM.Cmpd
    # Create a SAM plot of results
    mSet<-PlotSAM.Cmpd(mSet, "sam_imp_0_", 
                       format = "png", dpi=72, width=NA)
}