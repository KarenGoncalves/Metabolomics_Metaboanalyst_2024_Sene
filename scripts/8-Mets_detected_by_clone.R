# Mets present in each clone
library(tidyverse)
library(MetaboAnalystR)

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
inFiles <- paste0("Inputs/CleanUp_LCMSMS_", Analysis_modes, "_rawHeight.txt")
outFiles <- paste0("Results/Analytes_detected_by_clone_", Analysis_modes, ".txt")

for (i in 1:length(Analysis_modes)) {
    file = inFiles[i]
    output = out_files[i]
    
    
    input = read_delim(file)
    metadata = data.frame(Clone = input[1, -1] %>% unlist %>% unname ,
                          Replicate = names(input)[-1]) %>%
        filter(Clone != "QC")
    input_logical = as.data.frame(input[-1,],
                                  row.names = input[-1,1] %>% c)[-1] %>%
        apply(MARGIN=2, \(x) as.numeric(x) %>% as.logical )
    rownames(input_logical) = input$Sample[-1]
    
    groups = unique(metadata$Clone)
    
    sapply(groups, \(x) {
        rowSums(input_logical[, metadata$Replicate[metadata$Clone == x]]) > 2
    }) %>% as.data.frame %>%
        write_delim(file = output, delim = "\t", 
                    append = F, col_names = T)
}
