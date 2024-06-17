# Annotate deregulated analytes
library(tidyverse)
library(ggrepel)
#### Variables ####
Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
FDR_threshold = 0.05
FC_threshold = 1

#### Load data ####
load("Results/Siggenes_DAAs.RData")
annotation <- sapply(Analysis_modes, simplify = F, \(x) {
    paste0("Inputs/Corrected_LCMSMS_", x, "_identification.txt") %>%
        read_delim(delim = "\t", na = "null") %>%
        mutate(AnalysisMode = x)
}) %>% list_rbind() 

#### Keep only significantly DAAs ####
differential_abundance_sig <- 
    differential_abundance %>%
    list_rbind %>% filter(!is.na(pValue),
                          pValue < FDR_threshold,
                          abs(FoldChange) > FC_threshold
    ) %>%
    separate(col = Metabolite, 
             into = c("Rt", "Mz"), 
             sep = "/",  convert = T)

DAM_ids <- differential_abundance_sig %>%
    select(Rt, Mz, AnalysisMode) %>%
    unique

#### Get info on DAAs ####
annotation_DAAs <- 
    inner_join(annotation, DAM_ids,
               by = join_by("Average Rt(min)" == "Rt",
                            "Average Mz" == "Mz",
                            "AnalysisMode" == "AnalysisMode"))
Annotation_DAAs <- left_join(annotation_DAAs,
          differential_abundance_sig,
          by = join_by("Average Rt(min)" == "Rt",
                       "Average Mz" == "Mz",
                       "AnalysisMode" == "AnalysisMode")) %>%
    mutate(Clean_name = gsub("; (LC-|CE\\d).+$", "",
                             `Metabolite name`) %>%
               gsub(pattern="\\(*[Nn]ot validated.*",
                    replacement="") %>%
               gsub(pattern="^DDAO$",
                    replacement="Decyl dimethyl amine oxide") %>%
               gsub(pattern="^l-", replacement="L-") %>%
               gsub(pattern="^(NCGC\\d+-\\d+).+$",
                    replacement="\\1") %>%
               gsub(pattern="(2R,3S,4S,5R,6R)-6-(((4S,5aS,7S,11aR,12aS)-4,7-dihydroxy-3-((2R,5S)-5",
                    replacement="Cucurbitacin glycoside", fixed = T) %>%
               gsub(pattern="(2S,3S,4S,5S,6R)-2-(((2aR,4S,5aS,7S,11aR,12aS)-4-hydroxy-3",
                    replacement="Cyclocarposide", fixed = T) %>%
               gsub(pattern="(4aR,5aS,9R)-9-ethynyl-9a,11b-dimethylhexadecahydrocyclopenta[1,2]phenanthro",
                    replacement="Estrane steroid", fixed = T) %>%
               gsub(pattern="(R)-((2R,3S,4S,5R,6S)-6-((3-(2,3-dihydrobenzo[b][1,4]dioxin-6-yl)-4-oxo-4H-chromen-7-yl)oxy)-3,4,",
                    replacement="Isoflavonoid O-glycosides", fixed = T)
    )

# Replace with sentence case names
names_to_correct <- Annotation_DAAs$Clean_name %in% 
    c("THIAMINE", "ABIETIC ACID", "THIAMINE PYROPHOSPHATE")
Annotation_DAAs$Clean_name[names_to_correct] <- 
    str_to_sentence(Annotation_DAAs$Clean_name[names_to_correct])

# write_delim(Annotation_DAAs, 
#             "Results/differentially_abundant_analytes_annotation.txt",
#             delim = "\t", na = "NA")
# Annotation_DAAs %>% filter(`Metabolite name` != "Unknown") %>% 
#     select(`Metabolite name`, Clean_name, AnalysisMode, 
#            plotContrast, FoldChange, pValue, padj) %>% 
#     arrange(plotContrast, FoldChange) %>% View

# #### Plot proportion of annotated DAAs ####
# annotation_DAAs %>% 
#     mutate(Annotated = factor(ifelse(`Metabolite name` != "Unknown",
#                               "Yes", "No"),
#                               levels = c("Yes", "No")
#                               )
#            ) %>% 
#     group_by(Annotated, AnalysisMode) %>% 
#     summarize(Analytes = length(`Metabolite name`)) %>% 
#     group_by(AnalysisMode) %>%
#     summarize(Analytes = Analytes,
#               Annotated = Annotated,
#               Total = sum(Analytes)) %>% 
#     mutate(Proportion = Analytes / Total,
#            Mode = 
#                case_when(AnalysisMode == "HILIC_Positive" ~ "HILIC\nPositive",
#                          AnalysisMode == "RP_Positive" ~ "Reversed phase\nPositive",
#                          .default = "Reversed phase\nNegative")) %>%
#     ggplot(aes(x = "", Proportion, fill=Annotated,
#                label = Analytes)) + 
#     geom_col() + 
#     geom_text_repel(aes(color = Annotated),
#                     show.legend = F) +
#     facet_grid(~Mode, scales = "free_y") +
#     coord_polar(theta = "y", start = 0) +
#     scale_fill_manual(values = c("No" = "grey", "Yes" = "red")
#                       ) + 
#     scale_color_manual(values = c("No" = "grey30", "Yes" = "red")
#     ) + 
#     theme_void() +
#     theme(legend.position = "bottom")
# ggsave("plots/Proportion_annotated_DAAs.pdf",
#        height=2.5, width=5, dpi=1200)
# 
# 
#### Heatmap of annotated analytes ####
heatmap_data <- read_delim("Results/HeatmapData_longFormat.txt")  %>%
    separate(col = Metabolite, 
             into = c("Rt", "Mz"), 
             sep = "/",  convert = T) %>%
    inner_join(Annotation_DAAs[, c("Clean_name", "INCHIKEY",
                                   "AnalysisMode", "Formula",
                                   "Manually modified for annotation", 
                                   "Average Rt(min)", "Average Mz")] %>%
                   filter(Clean_name != "Unknown"),
               by = join_by("Rt" == "Average Rt(min)",
                            "Mz" == "Average Mz",
                            "AnalysisMode" == "AnalysisMode")) %>% 
    mutate(Metabolite_nameUnique = paste0(Clean_name, "\n", 
                                          Rt, "/", Mz)
    )

heatmap_data %>% dplyr::filter(Clean_name != "Unknown") %>% 
    dplyr::select(INCHIKEY, Clean_name) %>% 
    unique() %>%
    View()

order_mets = (heatmap_data %>%
                   dplyr::select(Contrast, FoldChange, Metabolite_nameUnique) %>%
                  unique() %>%
                   pivot_wider(names_from = Contrast, 
                               values_from = FoldChange,
                               values_fill=0) %>%
                   data.frame(row.names = .$Metabolite_nameUnique))[-1] %>%
                  as.matrix() %>% dist() %>% hclust

color_limit = max(abs(heatmap_data$FoldChange)) %>%
    ceiling()

color_scale = c(-color_limit, -color_limit/2,
                0, color_limit/2, color_limit)

heatmap_data %>%
    mutate(ordered_mets = Metabolite_nameUnique %>%
               factor(levels = order_mets$labels[order_mets$order]),
           Contrast_ordered = plotContrast %>%
               factor(levels = sapply(1:6, \(x) paste(tableContrast$Numerator[x], "vs",
                                                      tableContrast$Denominator[x]))
               )
    ) %>%
    ggplot(aes(y=ordered_mets, x = Contrast_ordered, fill = FoldChange)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), 
                         name = "Fold change",
                         limits = range(color_scale),
                         breaks = color_scale,
                         labels = color_scale) +
    labs(x = "", y="", fill="Fold Change") +
    theme_classic() +
    theme(legend.position = "bottom",
          #axis.text.x = element_blank(),
          #axis.ticks.x = element_blank()
          )
