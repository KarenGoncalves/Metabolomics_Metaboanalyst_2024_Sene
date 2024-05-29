# Annotate deregulated analytes
library(tidyverse)
library(MetaboAnalystR)


Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")

FC_res_long <- 
    read_delim("Results/FC_ANOVA_allmodes.txt", delim = "\t")

DAM_ids <- 
    FC_res_long %>% filter(Deregulation != "No") %>%
    select(Metabolite, AnalysisMode) %>%
    unique
