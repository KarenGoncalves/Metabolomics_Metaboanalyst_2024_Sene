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
               Denominator = c("pPTGE30", "pPTGE30", "pPTGE30", 
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
    differential_abundance_all %>%
    list_rbind() %>%
    filter(pValue < FDR_threshold,
           abs(FoldChange) > FC_threshold
    ) %>% 
    mutate(Regulation = case_when(FoldChange < -FC_threshold ~ "Down-regulated",
                                  FoldChange > FC_threshold ~ "Up-regulated",
                                  .default = "None")) 
write_delim(og_diff_abundance,
            "Results/HeatmapData_longFormat.txt",
            delim = "\t")

order_mets = ((og_diff_abundance %>%
                   dplyr::select(Contrast,FoldChange, Metabolite) %>%
                   pivot_wider(names_from = Contrast, values_from = FoldChange,
                               values_fill = 0) %>%
                   data.frame(row.names = .$Metabolite))[-1] %>%
                  as.matrix() %>% dist() %>% hclust)

color_limit = max(abs(og_diff_abundance$FoldChange)) %>%
    ceiling()

color_scale = c(-color_limit, -color_limit/2,
                0, color_limit/2, color_limit)
og_diff_abundance %>%
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

ggsave("plots/Siggenes_FC.pdf", 
    width=8, height=6)

ggsave("plots/Siggenes_FC.svg", 
       width=8, height=6)
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
    sapply(Analysis_modes, simplify = F, \(x) {
        (data = zscore_DAMs %>%
            filter(Analysis_mode == x) %>%
            dplyr::select(Metabolite, Replicate, zScore) %>%
            pivot_wider(id_cols = Metabolite,
                        names_from = Replicate,
                        values_from = zScore) %>%
            data.frame(row.names = .$Metabolite))[-1] %>%
    dist() %>% hclust()
    })
mets_per_AM = sapply(Analysis_modes, \(x) length(reorder_mets[[x]]$order))
Mets_ordered = data.frame(
    Metabolite = sapply(reorder_mets, \(x) x$labels) %>%
        unlist %>%
        unname(),
    order = c(reorder_mets[[1]]$order,
              reorder_mets[[2]]$order + mets_per_AM[1],
              reorder_mets[[3]]$order + mets_per_AM[1] + mets_per_AM[2]
    )
) %>% arrange(order)
nDAMs = nrow(Mets_ordered)

zscore_plots <- 
    sapply(Analysis_modes, simplify = F, \(x) {
        zscore_DAMs %>%
            filter(Analysis_mode == x) %>%
    mutate(ordered_Mets = Metabolite %>%
               factor(levels = reorder_mets[[x]]$labels[
                   reorder_mets[[x]]$order])
           )  %>%
    ggplot(aes(ordered_Mets, y = Replicate, fill = zScore)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), 
                         limits = c(-4, 4),
                         breaks = seq(-4, 4, 2),
                         labels = seq(-4, 4, 2)) +
    labs(x = "", y="", fill="Peak intensity\n(z-score)",
         title = gsub("_", " ", x)) +
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
})

library(ggpubr)
zscore_plots[2:3] <- sapply(
    zscore_plots[2:3], simplify = F, 
    \(x) {
        x + rremove('y.ticks') +
            rremove('y.axis') +
            rremove('y.text')
})
rel_widths = (mets_per_AM/sum(mets_per_AM)) + 
    c(0.1, 0, 0)
ggarrange(plotlist = zscore_plots, 
          vjust = T, common.legend = T,
          legend = "bottom", 
          ncol = 3, 
          widths = rel_widths[1]
)



ggsave("plots/Relative_abundance_DAAs.pdf", 
       height = 6, width = 8)
ggsave("plots/Relative_abundance_DAAs.svg", 
       height = 6, width = 8)
