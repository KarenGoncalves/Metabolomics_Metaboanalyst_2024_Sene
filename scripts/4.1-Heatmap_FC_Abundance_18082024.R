library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(RColorBrewer)
theme_classic() %>% theme_set()

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
inFiles <- paste0("Results/", Analysis_modes, "_mSet.RData")
mets_abundance = sapply(inFiles, simplify = F, \(x) {
    load(x)
    mSet$dataSet$filt %>% as.data.frame
})
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

#### Plot differential abundance ####

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


#### Plot peak normalized peak height ####
DAMs_abundance <- sapply(mets_abundance %>% names(), 
                         simplify = F, \(x) {
                             (t(mets_abundance[[x]]) %>% as.data.frame())[1:12] %>%
                                 as.data.frame() %>%
                                 mutate(Analysis_mode = gsub("_mSet.RData",
                                                             "", basename(x))
                                 )
                         }) %>% list_rbind() %>% 
    filter(rownames(.) %in% og_diff_abundance$Metabolite) %>%
    mutate(Metabolite = rownames(.)) %>%
    pivot_longer(cols = -c(Metabolite, Analysis_mode),
                 names_to="Replicate",
                 values_to="Height") %>%
    mutate(Clone = gsub("_I+$", "", Replicate)) 

summary_DAMs <- DAMs_abundance %>%
    group_by(Metabolite, Analysis_mode) %>%
    summarize(Mean_norm_abundance = mean(Height, na.rm=T),
              SD_norm_abundance = sd(Height, na.rm=T))

zscore_DAMs <- 
    sapply(1:nrow(DAMs_abundance), simplify = F, \(x) {
        Met=DAMs_abundance$Metabolite[x]
        info = summary_DAMs %>%
            filter(Metabolite == Met)
        DAMs_abundance[x,] %>%
            mutate(zScore = (Height - info$Mean_norm_abundance) /
                       info$SD_norm_abundance)
    }) %>% list_rbind()


reorder_mets <- 
        (zscore_DAMs %>%
             dplyr::select(Metabolite, Replicate, zScore) %>%
             pivot_wider(id_cols = Metabolite,
                         names_from = Replicate,
                         values_from = zScore) %>%
             data.frame(row.names = .$Metabolite))[-1,] %>%
            dist() %>% hclust()

zscore_plots <- 
        zscore_DAMs %>%
            mutate(ordered_Mets = Metabolite %>%
                       factor(levels = reorder_mets$labels[
                           reorder_mets$order])
            )  %>%
            ggplot(aes(ordered_Mets, y = Replicate, fill = zScore)) +
            geom_tile() +
            scale_fill_gradientn(colors = rev(brewer.pal(7, "RdBu")), 
                                 limits = c(-4, 4),
                                 breaks = seq(-4, 4, 2),
                                 labels = seq(-4, 4, 2)) +
            labs(x = "", y="", fill="Peak intensity\n(z-score)") +
            facet_grid(rows = vars(Clone), 
                       drop = T, scales = "free_y") +
            theme_classic() +
            theme(legend.position = "bottom",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  strip.background = element_blank(),
                  strip.text.y = element_blank(), 
                  panel.spacing = unit(-2, "mm"),
                  plot.margin = margin(r=-.2,
                                       t=0, b=0, l=-.3,
                                       'cm'))


ggsave(plot = zscore_plots, "plots/Relative_abundance_DAAs.pdf", 
       height = 6, width = 8)
