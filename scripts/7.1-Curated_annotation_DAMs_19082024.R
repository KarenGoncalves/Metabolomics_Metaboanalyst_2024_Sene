# Annotate deregulated analytes
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
#### Variables ####
Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
FDR_threshold = 0.05
FC_threshold = 1
empty_vector =  "pPTGE30"

#### Load data ####
load("Results/Siggenes_DAAs.RData")
annotation <- sapply(Analysis_modes, simplify = F, \(x) {
    paste0("Inputs/Corrected_LCMSMS_", x, "_identification.txt") %>%
        read_delim(delim = "\t", na = "null") %>%
        mutate(AnalysisMode = x)
}) %>% list_rbind() 

#### Keep only significantly DAAs ####
differential_abundance_sig <- 
    differential_abundance_all %>%
    list_rbind %>% filter(!is.na(pValue),
                          pValue < FDR_threshold,
                          abs(FoldChange) > FC_threshold
    ) %>% filter(grepl(empty_vector, Contrast)) %>%
    separate(col = Metabolite, 
             into = c("Rt", "Mz"), 
             sep = "/",  convert = T)

DAM_ids <- differential_abundance_sig %>%
    select(Rt, Mz, AnalysisMode) %>%
    unique
tableContrast = unique(differential_abundance_sig %>%
                           select(Contrast)) %>%
    separate(col = Contrast,
             into = c("Numerator", "Denominator"), sep = "_v_")

#### Get info on DAAs ####
Annotation_DAAs <- left_join(annotation,
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
               gsub(pattern="andrographolide",
                    replacement="Andrographolide", fixed=T) %>%
               gsub(pattern="(2R,3S,4S,5R,6R)-6-(((4S,5aS,7S,11aR,12aS)-4,7-dihydroxy-3-((2R,5S)-5-(2-hydroxypropan-2-yl)-2-methyltetrahydrofuran-2-yl)-2a,5a,8,8-tetramethylhexadecahydrocyclopenta[a]cyclopropa[e]phenanthren-9-yl)oxy)-5-(((2S,3R,4R)-3,4-dihydroxy-4-(hydroxymethyl)tetrahydrofuran-2-yl)oxy)-2-(hydroxymethyl)tetrahydro-2H-pyran-3,4-diol",
                    replacement="Cucurbitacin glycoside [2]", fixed=T) %>%
               gsub(pattern="(2S,3S,4S,5S,6R)-2-(((2aR,4S,5aS,7S,11aR,12aS)-4-hydroxy-3-((2R,5S)-5-(2-hydroxypropan-2-yl)-2-methyltetrahydrofuran-2-yl)-2a,5a,8,8-tetramethyl-9-(((2S,3R,4S,5R)-3,4,5-trihydroxytetrahydro-2H-pyran-2-yl)oxy)hexadecahydrocyclopenta[a]cyclopropa[e]phen...",
                    replacement="Cucurbitacin glycoside [1]", fixed=T) %>%
               gsub(pattern="(2S,3S,4S,5S,6R)-2-(((2aR,4S,5aS,7S,11aR,12aS)-4-hydroxy-3.+",
                    replacement="Cyclocarposide", fixed=T) %>%
               gsub(pattern="(4aR,5aS,9R)-9-ethynyl-9a,11b-dimethylhexadecahydrocyclopenta[1,2]phenanthro[8a,9-b]oxirene-3,9-diol",
                    replacement="Estrane steroid", fixed=T) %>%
               gsub(pattern="N-((octahydro-1H-quinolizin-1-yl)methyl)-2,4,5,6-tetrahydrocyclopenta[c]pyrazole-3-carboxamide",
                    replacement="Quinolizine", fixed=T) %>%
               gsub(pattern="2-(3,4-dimethoxyphenyl)-7-methoxy-4H-chromen-4-one",
                    replacement="7-O-methylated flavonoid [1]", fixed=T) %>%
               gsub(pattern="(R)-((2R,3S,4S,5R,6S)-6-((3-(2,3-dihydrobenzo[b][1,4]dioxin-6-yl)-4-oxo-4H-chromen-7-yl)oxy)-3,4,5-trihydroxytetrahydro-2H-pyran-2-yl)methyl 2-((tert-butoxycarbonyl)amino)-3-phenylpropanoate",
                    replacement="Isoflavonoid O-glycoside [1]", fixed=T) %>%
               gsub(pattern="methyl 2-((4-methyl-2-oxo-2H-chromen-7-yl)oxy)propanoate",
                    replacement="Coumarin derivative [1]", fixed=T) %>%
               gsub(pattern="4-((1R,3S,5r,7r)-5,7-dimethyl-1,3-diazaadamantan-2-yl)-2-methoxyphenol",
                    replacement="Methoxyphenol-type compound [1]", fixed=T) %>%
               gsub(pattern="7-benzyl-11,14-dimethyl-16-(2-methylpropyl)-10,13-di(propan-2-yl)-17-oxa-1,5,8,11,14-pentazabicyclo[17.3.0]docosane-2,6,9,12,15,18-hexone",
                    replacement="Cyclodepsipeptide", fixed=T) %>%
               gsub(pattern="7-hydroxy-3-(4-hydroxyphenyl)-8-((2S,3R,4R,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)tetrahydro-2H-pyran-2-yl)-4H-chromen-4-one",
                    replacement="Isoflavonoid C-glycoside [1]", fixed=T)           
           
    )

# Replace with sentence case names
names_to_correct <- Annotation_DAAs$Clean_name %in% 
    c("THIAMINE", "ABIETIC ACID", "THIAMINE PYROPHOSPHATE", "lappaconitine")
Annotation_DAAs$Clean_name[names_to_correct] <- 
    str_to_sentence(Annotation_DAAs$Clean_name[names_to_correct])

# write_delim(Annotation_DAAs, 
#             "Results/differentially_abundant_analytes_annotation.txt",
#             delim = "\t", na = "NA")
# Annotation_DAAs %>% filter(`Metabolite name` != "Unknown") %>% 
#     select(`Metabolite name`, Clean_name, AnalysisMode, 
#            plotContrast, FoldChange, pValue, padj) %>% 
#     arrange(plotContrast, FoldChange) %>% View


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
           Contrast_ordered = plotContrast %>%
               factor(levels = sapply(1:nrow(tableContrast), \(x) {
                   paste(tableContrast$Numerator[x], "vs",
                         tableContrast$Denominator[x])
               }))
    ) %>%
    ggplot(aes(y=ordered_mets, x = Contrast_ordered, fill = FoldChange)) +
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
          axis.text.x = element_text(angle=45, hjust=1, vjust=1),
          axis.text.y = element_text(size=7),
          # legend.key.height = unit(5, "mm"),
          # legend.key.width = unit(1, "cm"),
          # #axis.ticks.x = element_blank()
    )

ggsave("plots/Annotated_heatmap_allModes_19082024.pdf",
       height=8, width=7, dpi=1200)
