# MetaboAnalyst Analysis
library(tidyverse)
library(MetaboAnalystR)


in_file <- "Inputs/Combined_rawHeight.txt" 
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("pPTGE30", "pPTGE30", "pPTGE30", 
                               "AC9.3", "AC9.3", "AC9.2"))

mSet<-InitDataObjects(data.type = "pktable", 
                      anal.type = "stat", 
                      paired = FALSE);
mSet<-Read.TextData(mSetObj = mSet, 
                    filePath = in_file, 
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

### PCA ###
mSet$dataSet$meta.info = 
    mSet$dataSet$meta.info %>% 
    filter(Class != "QC")
mSet$dataSet$meta.info$Class <- 
    mSet$dataSet$meta.info$Class %>% 
    droplevels("QC")

mSet$dataSet$cls = 
    mSet$dataSet$meta.info$Class 
mSet$dataSet$norm = 
    mSet$dataSet$norm[1:nrow(mSet$dataSet$meta.info),]

#### Perform PCA analysis ####
plot_basename <- "plots/Combined_modes_"
mSet<-PCA.Anal(mSet)

mSet<-PlotPCAPairSummary(
    mSet, pc.num = 3, 
    imgName = paste0(plot_basename, "PCA_pair"),
    format = "png", dpi = 300, width=NA)
mSet$dataSet$table <- 
    mSet$dataSet$filt

mSet$dataSet$adjusted.mat <-
    mSet$dataSet$norm
mSet$dataSet$batch.cls <- mSet$dataSet$cls

# Create PCA overview
mSet<-PlotPCA.overview(
    mSet, 
    imgName = paste0(plot_basename, "PCA_overview"),
    format = "png", dpi = 300, method = "Pareto/Log")

# Create PCA scree plot
mSet<-PlotPCAScree(
    mSet, 
    imgName = paste0(plot_basename, "PCA_scree"),
    format = "png", dpi = 300, width=NA, scree.num = 5)

# Create a 2D PCA score plot
mSet<-PlotPCA2DScore(
    mSet,
    imgName = paste0(plot_basename, "PCA_score2d"), 
    format = "pdf", dpi=72, 
    pcx = 1, pcy = 2, reg = 0.99, 
    show = 0, grey.scale = 0)


# Create a 2D PCA score plot
mSet<-PlotPCABiplot(
    mSet,
    imgName = paste0(plot_basename, "PCA_biplot"), 
    format = "png", dpi=300, 
    width=NA, inx1 = 1, inx2 = 2)



