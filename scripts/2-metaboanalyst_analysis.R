# MetaboAnalyst Analysis
library(tidyverse)
library(MetaboAnalystR)

in_files = list.files(path = "Inputs/",
                      pattern = "CleanUp_",
                      full.names = T)
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("E30", "E30", "E30", 
                               "AC9.3", "AC9.3", "AC9.2"))

for (i in 1:3) {
    
    AnalysisMode <- gsub("CleanUp_LCMSMS_(.+)_rawHeight.txt",
                         "\\1", basename(in_files[i]))
    
    plots_basename <- paste0("plots/", AnalysisMode, "/")
    results_basename <- paste0("Results/", AnalysisMode)
    mSet<-InitDataObjects(data.type = "pktable", 
                          anal.type = "stat", 
                          paired = FALSE);
    mSet<-Read.TextData(mSetObj = mSet, 
                        filePath = in_files[i], 
                        format = "colu", 
                        lbl.type = "disc");
    mSet<-SanityCheckData(mSet);
    mSet<-ReplaceMin(mSet);
    mSet<-SanityCheckData(mSet);
    mSet<-FilterVariable(mSet, qc.filter = "T", rsd = 25, 
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
}

