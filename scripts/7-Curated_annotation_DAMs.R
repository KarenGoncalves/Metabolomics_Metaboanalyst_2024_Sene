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
left_join(annotation_DAAs,
          differential_abundance_sig,
          by = join_by("Average Rt(min)" == "Rt",
                       "Average Mz" == "Mz",
                       "AnalysisMode" == "AnalysisMode")) %>% 
    write_delim("Results/differentially_abundant_analytes_annotation.txt",
                delim = "\t", na = "NA")

#### Plot proportion of annotated DAAs ####
annotation_DAAs %>% 
    mutate(Annotated = factor(ifelse(`Metabolite name` != "Unknown",
                              "Yes", "No"),
                              levels = c("Yes", "No")
                              )
           ) %>% 
    group_by(Annotated, AnalysisMode) %>% 
    summarize(Analytes = length(`Metabolite name`)) %>% 
    group_by(AnalysisMode) %>%
    summarize(Analytes = Analytes,
              Annotated = Annotated,
              Total = sum(Analytes)) %>% 
    mutate(Proportion = Analytes / Total,
           Mode = 
               case_when(AnalysisMode == "HILIC_Positive" ~ "HILIC\nPositive",
                         AnalysisMode == "RP_Positive" ~ "Reversed phase\nPositive",
                         .default = "Reversed phase\nNegative")) %>%
    ggplot(aes(x = "", Proportion, fill=Annotated,
               label = Analytes)) + 
    geom_col() + 
    geom_text_repel(aes(color = Annotated),
                    show.legend = F) +
    facet_grid(~Mode, scales = "free_y") +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = c("No" = "grey", "Yes" = "red")
                      ) + 
    scale_color_manual(values = c("No" = "grey30", "Yes" = "red")
    ) + 
    theme_void() +
    theme(legend.position = "bottom")
ggsave("plots/Proportion_annotated_DAAs.pdf",
       height=2.5, width=5, dpi=1200)
