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
    mSet<-PCA.Anal(mSet)
    
    mSet<-PlotPCAPairSummary(
        mSet, pc.num = 3, 
        imgName = paste0(plot_basename[[i]], "PCA_pair_0_"),
        format = "png", dpi = 300, width=NA)
    
    # Create PCA scree plot
    mSet<-PlotPCAScree(
        mSet, 
        imgName = paste0(plot_basename[[i]], "PCA_scree_0_"),
        format = "png", dpi = 72, width=NA, scree.num = 5)
    
    # Create a 2D PCA score plot
    mSet<-PlotPCA2DScore(
        mSet,
        imgName = paste0(plot_basename[[i]], "PCA_score2d_0_"), 
        format = "png", dpi=72, 
        width=NA, pcx = 1, pcy = 2, reg = 0.95, 
        show = 0, grey.scale = 0)

    #### Partial least square - discriminant analysis (PLS-DA) ####
    
    mSet = runPLS.Analysis(
        mSetObj = mSet, 
        compNum = 8,
        iterations = 10000, 
        method = "oscorespls",  
        basenamePlots = plot_basename[[i]])
    # plot_RegCoef(mSet, 1)
    
    #### OPLS-DA ####
    mSet = OPLSR.Anal(mSetObj = mSet, reg = T)    
    mSet = OPLSDA.Permut(mSetObj = mSet, num = 1000)
    mSet = PlotOPLS.Permutation(mSet, 
                         imgName = paste0(plot_basename[[i]], 
                                          "OPLS_perm_2_"), 
                         format = "png", dpi=300, width=NA)
    # Create a 2D oPLS-DA score plot
    mSet<-PlotOPLS2DScore(mSet, 
                          imgName = paste0(plot_basename[[i]], 
                                           "OPLS_score2d_0_"), 
                          format = "png", 
                          dpi=72, 
                          width=NA, 
                          inx1 = 1, inx2 = 2,# PCs to plot
                          reg = 0.95, # conf.level of ellipses
                          show = 0, # show variable labels
                          grey.scale = 0)
    
    # Create a significant features plot
    mSet<-PlotOPLS.Splot(mSet, 
                         imgName = paste0(plot_basename[[i]], 
                                          "OPLS_splot_0_"), 
                         "all", "png", 300, width=NA);
    mSet$analSet$oplsda$splot.mat %>% 
        as.data.frame %>% 
        mutate(Feature = rownames(.)) %>% 
        ggplot(aes(x = `p[1]`, y = `p(corr)[1]`, label = Feature)) + 
        geom_point() + 
        geom_text_repel()
    ggsave(paste0(plot_basename[[i]], 
                  "OPLS_splot_0.png"),
           width = 6, height = 6, dpi = 300)
    # 
    # # Create a plot of the model overview
    mSet<-PlotOPLS.MDL(mSet, 
                       imgName = paste0(plot_basename[[i]], 
                                        "OPLS_mdl_0_"), 
                       format = "png", 
                       dpi=72, width=NA)
    save(mSet, file = inFiles[i])
}
