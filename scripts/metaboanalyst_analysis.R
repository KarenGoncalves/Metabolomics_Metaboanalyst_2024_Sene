# MetaboAnalyst Analysis
# install.packages("agricolae")
library(tidyverse)
library(MetaboAnalystR)
library(car)
library(dunn.test)
library(agricolae)
source("scripts/help_script/ANOVA_mets.R")

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
    plots_basename <- paste0("plots/", AnalysisMode)
    results_basename <- paste0("Results/", AnalysisMode)
    mSet<-InitDataObjects("conc", "stat", FALSE);
    mSet<-Read.TextData(mSet, in_files[i], "colu", "disc");
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
                                                 "_rowNormalization_"), 
                          format ="png", dpi=300, width=NA);
    mSet<-PlotSampleNormSummary(mSet, imgName = paste0(plots_basename,
                                                       "_colNormalization_"), 
                                format ="png", dpi=300, width=NA);
    
    # Perform fold-change analysis on uploaded data, unpaired
    mSet$analSet$fc<- MetaboAnalystR::GetFC(mSetObj = mSet, paired = F, 
                                            tableContrast = tableContrast)
    
    # Plot fold-change analysis
    FC_res_table <- data.frame(IDs = mSet$analSet$fc$AC9.1_v_E30$fc.log %>% 
                                   names)
    for (i in names(mSet$analSet$fc)) {
        FC_res_table[[i]] = mSet$analSet$fc[[i]]$fc.log
    }
    
    FC_res_long <- FC_res_table %>% pivot_longer(cols = !IDs,
                                                 values_to = "log2FC",
                                                 names_to = "Contrasts") %>%
        mutate(plotContrast = gsub("_v_", " vs ", Contrasts),
               Regulation = case_when(log2FC > 2 ~ "Up-regulated",
                                      log2FC < -2 ~ "Down-regulated",
                                      .default = "No")) %>%
        separate(col = Contrasts, into = c("Clone1", "Clone2"),
                 sep = "_v_") %>%
        rename("Metabolite" = IDs)
    
    
    FC_res_long %>% ggplot(aes(y = log2FC, x = Metabolite, 
                               color = Regulation, alpha = Regulation)) +
        geom_point() + 
        scale_color_manual(values = c("Up-regulated" = "red",
                                      "Down-regulated" = "blue",
                                      "No" = "grey50"),
                           name = "") +
        scale_alpha_manual(values = c("Up-regulated" = 1,
                                      "Down-regulated" = 1,
                                      "No" = 0.3), 
                           name = "") +
        facet_wrap(~plotContrast, nrow = 3, ncol = 2, as.table = F) +
        theme_classic() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom") +
        labs(x = "Analytes", y = bquote("log"[2]~"FC"))
    ggsave(paste0(plots_basename, "_FCAnalysis.png"), dpi = 300,
           height = 7, width = 6)
    
    # To view fold-change 
    mets <- names(mSet$dataSet$norm)
    clones <- mSet$dataSet$meta.info %>%
        filter(Class != "QC")
    data <- mSet$dataSet$norm[1:nrow(clones), ] %>% as.data.frame()
    anova_result <- ANOVA_mets(data, clones, mets)
    
    heteroscedastic_mets <- sapply(anova_result, \(x) x$homoscedascity < .05) %>%
        which %>% length
    DAM_tables <- sapply(mets, simplify = F, \(x) {
        anova_result[[x]]$posthoc
    }) %>% list_rbind
    
    write_delim(DAM_tables, 
                paste0(results_basename, "_DAM_table.txt"),
                delim = "\t")
    write_delim(FC_res_long, 
                paste0(results_basename, "_FC_result.txt"))
    save(mSet, file = paste0(results_basename, "_mSet.RData"))
}