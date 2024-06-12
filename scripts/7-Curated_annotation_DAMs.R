# Annotate deregulated analytes
library(tidyverse)
library(MetaboAnalystR)

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

# write_delim(annotation_DAAs,
#             file="Results/annotation_DAAs.txt",
#             delim="\t", na="NA")
