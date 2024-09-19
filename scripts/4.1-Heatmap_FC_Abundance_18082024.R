library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(RColorBrewer)
theme_classic() %>% theme_set()

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
clones = c(paste0("AC9.", 1:3), "pPTGE30")
plot_basename <- paste0("plots/", Analysis_modes, "/")
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("pPTGE30", "pPTGE30", "pPTGE30", 
                               "AC9.3", "AC9.3", "AC9.2"))
FDR_threshold = 0.05
FC_threshold = 1

#### Load data ####
kept_mets = sapply(Analysis_modes, simplify = F, \(x) {
    load(paste0("Results/", x, "_mSet.RData"))
    data.frame(MetID=mSet$dataSet$filt %>% colnames,
               AnalysisMode = x)
}) %>% list_rbind() %>%
    separate(MetID, into = c("Rt", "Mz"), convert = T, sep="/")

mets_abundance <- 
    sapply(Analysis_modes, simplify = F, \(x) {
        read_delim(paste0("Inputs/CleanUp_LCMSMS_", x, "_rawHeight.txt"))[-1,] %>%
            mutate(AnalysisMode = x) %>%
            separate(Sample, into = c("Rt", "Mz"), convert = T, sep="/") 
}) %>% list_rbind

mets_abundance <- mets_abundance[, grep("QC", names(mets_abundance), invert = T)] %>%
    inner_join(kept_mets,
               by = c("Rt", "Mz", "AnalysisMode")
    )
cols.num <- sapply(clones, \(x) grep(x, names(mets_abundance), value = T)) %>%
    c
mets_abundance[cols.num] <- sapply(mets_abundance[cols.num],as.numeric)
sapply(mets_abundance, class)


#### Plot differential abundance ####
load("Results/Siggenes_DAAs.RData")

og_diff_abundance <- 
    differential_abundance_all %>%
    list_rbind() %>%
    filter(pValue < FDR_threshold,
           abs(FoldChange) > FC_threshold
    ) %>% 
    mutate(Regulation = case_when(FoldChange < -FC_threshold ~ "Down-regulated",
                                  FoldChange > FC_threshold ~ "Up-regulated",
                                  .default = "None")) 


diff_abundance_pPTGE30 <- 
    og_diff_abundance %>%
    filter(grepl("pPTGE30", Contrast))


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

order_mets = ((diff_abundance_pPTGE30 %>%
                   dplyr::select(Contrast,FoldChange, Metabolite) %>%
                   pivot_wider(names_from = Contrast, values_from = FoldChange,
                               values_fill = 0) %>%
                   data.frame(row.names = .$Metabolite))[-1] %>%
                  as.matrix() %>% dist() %>% hclust)

color_limit = max(abs(diff_abundance_pPTGE30$FoldChange)) %>%
    ceiling()

color_scale = c(-color_limit, -color_limit/2,
                0, color_limit/2, color_limit)
diff_abundance_pPTGE30 %>%
    mutate(ordered_mets = Metabolite %>%
               factor(levels = order_mets$labels[order_mets$order]),
           Contrast_ordered = plotContrast %>%
               factor(levels = sapply(1:3, \(x) paste(tableContrast$Numerator[x], "vs",
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

ggsave("plots/Siggenes_FC_pPTGE30.pdf", 
       width=8, height=6)

#### Plot peak normalized peak height ####
DAMs_abundance <- mets_abundance %>%
    rename("Analysis_mode" = AnalysisMode) %>%
    mutate(Metabolite = paste0(Rt, "/", Mz)) %>%
    filter(Metabolite %in% og_diff_abundance$Metabolite) %>%
    pivot_longer(cols = -c(Metabolite, Analysis_mode, Rt, Mz),
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
            mutate(
                zScore = (Height - info$Mean_norm_abundance) / info$SD_norm_abundance,
                relAbundance = Height/mean(info$Mean_norm_abundance) %>% log2)
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


ggsave(plot = zscore_plots, paste0("plots/Relative_abundance_DAAs_", Sys.Date(), ".pdf"), 
       height = 6, width = 8)

#### All up/down AC9
Metabolites_AC9 <- 
    (og_diff_abundance %>%
         filter(grepl("pPTGE30", Contrast)) %>%
         group_by(Metabolite) %>%
         summarise(`Down-regulated` = length(which(Regulation == "Down-regulated")),
                   `Up-regulated` = length(which(Regulation == "Up-regulated"))) %>%
         filter(`Up-regulated` == 3 | `Down-regulated` == 3))[[1]]

order_mets = ((diff_abundance_pPTGE30 %>%
                   filter(Metabolite %in% Metabolites_AC9) %>% 
                   dplyr::select(Contrast,FoldChange, Metabolite) %>%
                   pivot_wider(names_from = Contrast, values_from = FoldChange,
                               values_fill = 0) %>%
                   data.frame(row.names = .$Metabolite))[-1] %>%
                  as.matrix() %>% dist() %>% hclust)

color_limit = max(abs(og_diff_abundance$FoldChange)) %>%
    ceiling()

color_scale = c(-color_limit, -color_limit/2,
                0, color_limit/2, color_limit)

    
diff_abundance_pPTGE30 %>%
    filter(Metabolite %in% Metabolites_AC9) %>%
    mutate(ordered_mets = Metabolite %>%
               factor(levels = order_mets$labels[order_mets$order]),
           Contrast_ordered = plotContrast %>%
               factor(levels = sapply(1:3, \(x) paste(tableContrast$Numerator[x], "vs",
                                                      tableContrast$Denominator[x]))
               )
    ) %>%
    ggplot(aes(ordered_mets, x = Contrast_ordered, 
               fill = FoldChange)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), 
                         limits = range(color_scale),
                         breaks = color_scale,
                         labels = color_scale) +
    labs(x = "", y="", 
         fill="Fold\nChange") +
    theme_classic() +
    theme(axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7),
          axis.ticks.x = element_blank())

ggsave("plots/FoldChange_AllUp_or_AllDown_AC9.pdf", width=5, height = 9,
    dpi=1200)

pdf("plots/AllUp_or_AllDown_AC9.pdf", width=8, height = 8)
zscore_DAMs %>%
    filter(Metabolite %in% Metabolites_AC9)  %>%
    ggplot(aes(y=Metabolite, x = Replicate, fill = zScore)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(7, "RdBu")), 
                         limits = c(-4, 4),
                         breaks = seq(-4, 4, 2),
                         labels = seq(-4, 4, 2)) +
    labs(x = "", y="", fill="Relative peak intensity\n(z-score)",
         caption = "zscore = x - mean / std. dev") +
    facet_grid(cols = vars(Clone), 
               drop = T, scales = "free_x") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.text.y=element_text(size=7),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(-2, "mm"))

pctDAMs <- 
    sapply(DAMs_abundance$Metabolite %>% unique, 
           simplify = F, \(Met) {
               data = DAMs_abundance %>%
                   filter(Metabolite == Met)
               maxAbundance = max(data$Height)
               data %>%
                   mutate(relAbundance = Height * 100 / maxAbundance)
           }) %>% list_rbind()

pctDAMs %>%
    filter(Metabolite %in% Metabolites_AC9)  %>% 
ggplot(aes(y=Metabolite, x = Replicate, fill = relAbundance)) +
    geom_tile() +
    scale_fill_gradient(low = "white", 
                        high = "#99000D") +
    labs(x = "", y="", fill="Relative peak intensity",
         caption = "rel peak intensity = log2(x * mean)") +
    facet_grid(cols = vars(Clone), 
               drop = T, scales = "free_x") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.text.y=element_text(size=7),
          axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(-2, "mm"))
dev.off()

