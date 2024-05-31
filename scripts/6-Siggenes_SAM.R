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
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("E30", "E30", "E30", 
                               "AC9.3", "AC9.3", "AC9.2"))
FDR_threshold = 0.001
adjusted_FDR_threshold = 0.000165

differential_abundance <- list()
for (i in 1:3) {
    pdf(file = paste0(plot_basename[i], 
                      "Differential_analysis.pdf"),
        width = 6, height = 6, onefile = T)
    rm("mSet")
    load(inFiles[i])
    differential_abundance[[i]] <- list()
    for (j in 1:nrow(tableContrast)){
        print(tableContrast[j,])
        numerator = tableContrast$Numerator[j]
        denominator = tableContrast$Denominator[j]
        contrastName = paste0(numerator, "_v_", denominator)
        replicates = rownames(mSet$dataSet$norm)
        mSet_subset_norm = mSet$dataSet$norm[c(grep(numerator, replicates),
                                               grep(denominator, replicates)), ]
        cls_subset = factor(rep(tableContrast[j,] %>% unlist, each = 3),
                            levels = tableContrast[j,])
        
        SAM = SAM = siggenes::sam(t(mSet_subset_norm), 
                                  cl = cls_subset, 
                                  method = "d.stat", B=100, 
                                  gene.names = names(mSet_subset_norm), 
                                  med = F, use.dm = T, 
                                  var.equal = T,
                                  control = samControl(n.delta = 10),
                                  R.fold = 2)
        SAM@msg = c(contrastName, SAM@msg)
        
        deltaValue = findDelta(SAM, fdr = adjusted_FDR_threshold)
        if (nrow(deltaValue) %>% is.null) {
            deltaValue = deltaValue["Delta"]
        } else {
            deltaValue = deltaValue[order(deltaValue[,"FDR"]), "Delta"]
        }
        sam.plot2(SAM, deltaValue, sig.col = c("blue", "red"), 
                  main = paste0(numerator, " vs ", denominator)
        )
        
        differential_abundance[[i]][[contrastName]] <- 
            data.frame(Analysis_type = Analysis_modes[i],
                       Contrast = contrastName,
                       plotContrast = paste0(numerator, " vs ", denominator),
                       Metabolites = names(SAM@d),
                       FC = SAM@fold,
                       log2FC = log2(SAM@fold),
                       pValue = SAM@p.value) %>%
            pivot_longer(cols = c(FC, log2FC, pValue),
                         names_to = "ResultType",
                         values_to = "Value")
                         
    }
    differential_abundance[[i]] <- 
        differential_abundance[[i]] %>%
        list_rbind() %>%
        pivot_wider(names_from = ResultType,
                    values_from = Value)
    # mSet <- SAM.Anal(mSetObj = mSet, 
    #                  method = "d.stat", # d.stat for parametric 
    #                  paired = F,
    #                  varequal = T, 
    #                  imgName = paste0(plot_basename[i], 
    #                                   "Differential_analysis_"),
    #                  dpi=300)
    # mSet<-PlotSAM.FDR(mSet, 
    #                   imgName = paste0(plot_basename[i], 
    #                                    "Differential_analysis_FDR"),
    #                   format = "png", 
    #                   dpi=300, width=NA)
    # ?PlotSAM.Cmpd
    # # Create a SAM plot of results
    # mSet<-PlotSAM.Cmpd(mSet, 
    #                    imgName = paste0(plot_basename[i], 
    #                                     "Differential_analysis_ImportantFeaures_"),
    #                    format = "png", dpi=300, width=NA)
dev.off()
}
dev.off()

subset_diff_abundance <- 
    differential_abundance %>%
    list_rbind() %>%
    filter(!is.na(FC), pValue < .05,
           #!grepl("AC9.2", Contrast)
           )

order_mets = ((subset_diff_abundance %>%
    dplyr::select(Contrast,log2FC, Metabolites) %>%
    pivot_wider(names_from = Contrast, values_from = log2FC,
                values_fill = 0) %>%
    data.frame(row.names = .$Metabolites))[-1] %>%
    as.matrix() %>% dist() %>% hclust)

subset_diff_abundance %>%
    mutate(ordered_mets = Metabolites %>%
               factor(levels = order_mets$labels[order_mets$order]),
           Contrast_ordered = plotContrast %>%
               factor(levels = sapply(1:6, \(x) paste(tableContrast$Numerator[x], "vs",
                                                      tableContrast$Denominator[x]))
                      )
           ) %>%
    ggplot(aes(ordered_mets, y = Contrast_ordered, fill = log2FC)) +
    geom_tile() +
    scale_fill_gradient2(high = "red",
                         mid="white",
                         low="blue",
                         na.value = "white") +
    labs(x = "", y="", fill="Fold Change") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
                          
