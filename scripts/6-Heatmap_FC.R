# Heatmaps
library(tidyverse)
library(MetaboAnalystR)
library(RColorBrewer)


Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("E30", "E30", "E30", 
                               "AC9.3", "AC9.3", "AC9.2"))

FC_res_long <- 
    read_delim("Results/FC_ANOVA_allmodes.txt", delim = "\t")

FC_res_long$plotContrast <- 
    FC_res_long$plotContrast %>% 
    factor(levels = 
               apply(tableContrast, 1, \(x) paste(x[1], x[2], 
                                                  sep = " vs "))
    )

Mets_heatmap_data <- 
    FC_res_long %>% filter(Deregulation != "No") 
     
metabolites_dendrogram <- 
    (Mets_heatmap_data %>%
    select(plotContrast, Metabolite, log2FC) %>%
    pivot_wider(id_cols = Metabolite, 
                names_from = plotContrast,
                values_from = log2FC, values_fill = 0) %>%
    data.frame(row.names = .$Metabolite))[-1] %>%
    dist() %>% hclust

Mets_heatmap_data$Metabolite <- 
    Mets_heatmap_data$Metabolite %>% 
    factor(levels = metabolites_dendrogram$labels[
        metabolites_dendrogram$order
    ])

max(abs(Mets_heatmap_data$log2FC))
Mets_heatmap_data %>%
    ggplot(aes(x=Metabolite, 
               y = plotContrast,
               fill = log2FC)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), 
                         name = "Fold change",
                         limits = c(-18, 18),
                         breaks = c(-18, -9, 0, 9, 18),
                         labels = c(-18, -9, 0, 9, 18)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          axis.ticks = element_blank())
ggsave("plots/Heatmap_metabolites_FC2.png",
       dpi=300, width=8, height=5)
    