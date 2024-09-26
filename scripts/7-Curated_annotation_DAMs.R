# Annotate deregulated analytes
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
#### Variables ####
Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
FDR_threshold = 0.01
FC_threshold = 1
empty_vector =  "pPTGE30"

#### Load data ####
annotation <- read_delim("Results/annotation_DAAs.txt")


tableContrast = unique(annotation %>%
                           select(Contrast)) %>%
    separate(col = Contrast,
             into = c("Numerator", "Denominator"), sep = "_v_")

annotation %>%
    select(Contrast, FoldChange, Clean_name) %>%
    filter(Clean_name != "Unknown") %>%
    unique()


#### Heatmap of annotated analytes ####
heatmap_data <- 
    annotation %>%
    rename("Rt" = "Average Rt(min)",
           "Mz" = "Average Mz",
           "AnalysisMode" = "AnalysisMode") %>% 
    mutate(Metabolite_nameUnique = paste0(Clean_name, " ", 
                                          Rt, "/", Mz)
    ) %>% filter(grepl(empty_vector, Contrast))  %>% 
    filter(Clean_name != "Unknown" &
               !grepl("w/o MS2", Clean_name)) 

putative_annotations <- 
    heatmap_data %>% select(Clean_name, INCHIKEY, Metabolite_nameUnique) %>%
    unique %>% group_by(Clean_name) %>% 
    summarize(uniqueName = Metabolite_nameUnique,
              n = length(Metabolite_nameUnique)) %>%
    mutate(namePlot = 
               ifelse(n > 1, 
                      paste0("Putative: ", uniqueName),
                      uniqueName))

heatmap_data$Metabolite_nameplot <- 
    sapply(heatmap_data$Metabolite_nameUnique,
           \(x) {
               putative_annotations$namePlot[
                   putative_annotations$uniqueName == x
               ]
           })

order_mets = (heatmap_data %>%
                  dplyr::select(Contrast, FoldChange, Metabolite_nameplot) %>%
                  unique() %>%
                  pivot_wider(names_from = Contrast, 
                              values_from = FoldChange,
                              values_fill=0) %>%
                  data.frame(row.names = .$Metabolite_nameplot))[-1] %>%
    as.matrix() %>% dist() %>% hclust

color_limit = max(abs(heatmap_data$FoldChange)) %>%
    ceiling()
color_scale = c(-color_limit, -color_limit/2, #-FC_threshold,
                0, #FC_threshold, 
                color_limit/2, color_limit)
colors = rev(brewer.pal(5, "RdBu"))
names(colors) = color_scale
colors["0"] = "white"

heatmap_data %>%
    mutate(ordered_mets = Metabolite_nameplot %>%
               factor(levels = order_mets$labels[order_mets$order]),
           strainsInterest = gsub(" .+", "", plotContrast)
    ) %>%
    ggplot(aes(y=ordered_mets, x = strainsInterest, fill = FoldChange)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors, 
                         name = "Fold change",
                         limits = range(color_scale),
                         breaks = color_scale,
                         labels = color_scale) +
    labs(x = "", y="", fill="Fold Change") +
    theme_classic() +
    theme(
        # legend.position = "bottom",
          axis.text = element_text(size=10),
          # legend.key.height = unit(5, "mm"),
          # legend.key.width = unit(1, "cm"),
          # #axis.ticks.x = element_blank()
    )

ggsave("plots/Annotated_heatmap_allModes_26092024.pdf",
       height=8, width=7, dpi=1200)
