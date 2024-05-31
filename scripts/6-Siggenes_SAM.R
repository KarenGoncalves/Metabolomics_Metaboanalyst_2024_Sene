# MetaboAnalyst multivariate analyses
# install.packages(c("ellipse", "pls"))
library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(RColorBrewer)
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
met_abundance <- list()
differential_abundance <- list()
for (i in 1:3) {
    pdf(file = paste0(plot_basename[i], 
                      "Differential_analysis.pdf"),
        width = 6, height = 6, onefile = T)
    rm("mSet")
    load(inFiles[i])
    met_abundance[[i]] <- mSet$dataSet$norm 
    differential_abundance[[i]] <- list()
    
    for (j in 1:nrow(tableContrast)){
        print(tableContrast[j,])
        numerator = tableContrast$Numerator[j]
        denominator = tableContrast$Denominator[j]
        contrastName = paste0(numerator, "_v_", denominator)
        replicates = rownames(mSet$dataSet$norm)
        mSet_subset_norm = met_abundance[[i]][c(grep(numerator, replicates),
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

color_limit = max(abs(subset_diff_abundance$log2FC)) %>%
    ceiling()

color_scale = c(-color_limit, -color_limit/2,
                0, color_limit/2, color_limit)
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
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), 
                         name = "Fold change",
                         limits = range(color_scale),
                         breaks = color_scale,
                         labels = color_scale) +
    labs(x = "", y="", fill="Fold Change") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
                          

DAMs_abundance <- sapply(met_abundance, simplify = F, \(x) {
    (t(x) %>% as.data.frame())[1:12]
}) %>% list_rbind() %>% 
    filter(rownames(.) %in% order_mets$labels) %>%
    mutate(Metabolite = factor(rownames(.),
                               levels = order_mets$labels)) %>%
    pivot_longer(cols = !Metabolite,
                 names_to="Replicate",
                 values_to="Norm_abundance") %>%
    mutate(Clone = gsub("_I+$", "", Replicate)) %>%
    filter(Metabolite %in% order_mets$labels) 

summary_DAMs <- DAMs_abundance %>%
    group_by(Metabolite) %>%
    summarize(Mean_norm_abundance = mean(Norm_abundance, na.rm=T),
              SD_norm_abundance = sd(Norm_abundance, na.rm=T))

zscore_DAMs <- 
    sapply(1:nrow(DAMs_abundance), simplify = F, \(x) {
        Met=DAMs_abundance$Metabolite[x]
        info = summary_DAMs %>%
            filter(Metabolite == Met)
        DAMs_abundance[x,] %>%
            mutate(zScore = (Norm_abundance - info$Mean_norm_abundance) /
                       info$SD_norm_abundance)
    }) %>% list_rbind()


reorder_mets <- (zscore_DAMs %>% 
    dplyr::select(Metabolite, Replicate, zScore) %>%
    pivot_wider(id_cols = Metabolite,
                names_from = Replicate,
                values_from = zScore) %>%
    data.frame(row.names = .$Metabolite))[-1] %>%
    dist() %>% hclust()


nDAMs = length(reorder_mets$order)
zscore_DAMs %>%
    mutate(ordered_Mets = Metabolite %>%
               factor(levels = reorder_mets$labels[reorder_mets$order])) %>%
    ggplot(aes(ordered_Mets, y = Replicate, fill = Norm_abundance)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), 
                         limits = c(-4, 4),
                         breaks = seq(-4, 4, 2),
                         labels = seq(-4, 4, 2)) +
    labs(x = "", y="", fill="Normalized\nabundance",
         caption=paste("Normalized abundance of deregulated analytes (all three analysis modes).
                       Deregulated analytes:", nDAMs)) +
    facet_wrap(~Clone, nrow = 4, drop = T, scales = "free_y") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(), 
          panel.spacing.y = unit(-2, "mm"))

ggsave("plots/Relative_abundance_DAAs.pdf", 
       height = 6, width = 8)
