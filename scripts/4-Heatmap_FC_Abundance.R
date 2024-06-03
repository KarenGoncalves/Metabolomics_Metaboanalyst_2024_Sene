library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(RColorBrewer)
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
FDR_threshold = 0.05
FC_threshold = 1
# FoldChanges <- read_delim(paste0("Results/FC_ANOVA_allmodes.txt")) %>%
#     dplyr::select(Metabolite,
#                   Clone1, Clone2,
#                   log2FC, plotContrast,
#                   AnalysisMode)
load("Results/Siggenes_DAAs.RData")

og_diff_abundance <- 
    differential_abundance %>%
    list_rbind() %>%
    filter(pValue < FDR_threshold,
           abs(FoldChange) > FC_threshold
    ) %>% 
    mutate(Regulation = case_when(FoldChange < -FC_threshold ~ "Down-regulated",
                                  FoldChange > FC_threshold ~ "Up-regulated",
                                  .default = "None")) 

subsets_diff_abundance <- 
    sapply(paste0("AC9.", 1:3), simplify = F, \(x) {
        og_diff_abundance %>%
            filter(!grepl(x, plotContrast))
    })

pdf("plots/Siggenes_FC.pdf", onefile = T,
    width=8, height=6)

list("Original" = og_diff_abundance,
     "AC9.2_AC9.3" = subsets_diff_abundance[[1]],
     "AC9.1_AC9.3" = subsets_diff_abundance[[2]],
     "AC9.1_AC9.2" = subsets_diff_abundance[[3]]) %>% 
    sapply(simplify = F, \(x) {
        order_mets = ((x %>%
                           dplyr::select(Contrast,FoldChange, Metabolite) %>%
                           pivot_wider(names_from = Contrast, values_from = FoldChange,
                                       values_fill = 0) %>%
                           data.frame(row.names = .$Metabolite))[-1] %>%
                          as.matrix() %>% dist() %>% hclust)
        
        color_limit = max(abs(x$FoldChange)) %>%
            ceiling()
        
        color_scale = c(-color_limit, -color_limit/2,
                        0, color_limit/2, color_limit)
        x %>%
            mutate(ordered_mets = Metabolite %>%
                       factor(levels = order_mets$labels[order_mets$order]),
                   Contrast_ordered = plotContrast %>%
                       factor(levels = sapply(1:6, \(x) paste(tableContrast$Numerator[x], "vs",
                                                              tableContrast$Denominator[x]))
                       )
            ) %>%
            ggplot(aes(ordered_mets, y = Contrast_ordered, fill = FoldChange)) +
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
        
    })
dev.off()

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
