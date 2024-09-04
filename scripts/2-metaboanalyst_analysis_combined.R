# MetaboAnalyst Analysis
library(tidyverse)
library(MetaboAnalystR)

in_files = "Inputs/Combined_rawHeight.txt"
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("pPTGE30", "pPTGE30", "pPTGE30", 
                               "AC9.3", "AC9.3", "AC9.2"))

plots_basename <- "plots/"
results_basename <- "Results/"
mSet<-InitDataObjects(data.type = "pktable", 
                      anal.type = "stat", 
                      paired = FALSE);
mSet<-Read.TextData(mSetObj = mSet, 
                    filePath = in_files, 
                    format = "colu", 
                    lbl.type = "disc");
mSet<-SanityCheckData(mSet);
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet);
mSet<-FilterVariable(mSet, qc.filter = "T", rsd = 10, 
                     var.filter = "none", var.cutoff = -1, 
                     int.filter = "mean", int.cutoff = 0)
mSet<-PreparePrenormData(mSet);
mSet<-Normalization(mSet, "GroupPQN", transNorm = "LogNorm", 
                    scaleNorm = "ParetoNorm", ref = "QC", 
                    ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, imgName = paste0(plots_basename,
                                             "rowNormalization_"), 
                      format ="png", dpi=300, width=NA);
mSet<-PlotSampleNormSummary(mSet, imgName = paste0(plots_basename,
                                                   "colNormalization_"), 
                            format ="png", dpi=300, width=NA);
save(mSet, file = paste0(results_basename, "_mSet.RData"))

