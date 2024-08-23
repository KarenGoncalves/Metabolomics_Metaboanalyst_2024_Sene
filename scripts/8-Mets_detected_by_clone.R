# Mets present in each clone
library(tidyverse)
library(MetaboAnalystR)
library(venn)

clones <- c(paste0("AC9.", 1:3), "pPTGE30")
Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
inFiles <- paste0("Inputs/CleanUp_LCMSMS_", Analysis_modes, "_rawHeight.txt")
outFiles <- paste0("Results/Analytes_detected_by_clone_", Analysis_modes, ".txt")
mSetData <- paste0("Results/", Analysis_modes, "_mSet.RData")

analytes_detected <- sapply(1:length(Analysis_modes), simplify = F, \(i) {
    file = inFiles[i]
    output = outFiles[i]
    load(mSetData[i])
    
    featuresKept_metaboanalyst <- 
        mSet[["dataSet"]][["filt"]] %>% colnames
    
    input = read_delim(file)
    names(input) <- gsub("^E30", "pPTGE30", names(input))
    input[1, -1] <- data.frame(gsub("^E30", "pPTGE30", input[1, -1])) %>% t
    metadata = data.frame(Clone = input[1, -1] %>% unlist %>% unname ,
                          Replicate = names(input)[-1]) %>%
        filter(Clone != "QC")
    input_logical = as.data.frame(input[-1,],
                                  row.names = input[-1,1] %>% c)[-1] %>%
        apply(MARGIN=2, \(x) as.numeric(x) %>% as.logical )
    rownames(input_logical) = input$Sample[-1]
    
    groups = unique(metadata$Clone)
    
    sapply(groups, \(x) {
        rowSums(input_logical[featuresKept_metaboanalyst,
                              metadata$Replicate[metadata$Clone == x]]) > 2
    }) %>% as.data.frame
}) %>% list_rbind() %>%
    filter(rowSums(.) != 0) %>%
    mutate(MetID = row.names(.)) %>%
    separate(col = MetID, into = c("RT", "Mz"), sep = "/", convert = T)

Annotation <- 
    read_delim("Results/Annotation_clean_names.txt") %>%
    right_join(analytes_detected, 
               by = join_by("Average Rt(min)" == "RT",
                            "Average Mz" == "Mz")) 

annotated_detected <- Annotation %>%
    filter(Clean_name != "Unknown" &
               !grepl("w/o MS2", Clean_name))


#### Venn plots ####
par(cex=1, font=1); 
venn_detected <- analytes_detected %>%
    select(all_of(clones))

svg('plots/Venn_filteredMetaboanalyst_detected.svg',
    width=8, height=8, pointsize = 30)
venn(venn_detected, 
     ilabels = "counts", 
     ellipse = T, box = F, 
     zcolor = "style", 
     par = T)
dev.off()

venn_annotated_detected <- 
    annotated_detected %>%
    select(all_of(clones))


svg('plots/Venn_filteredMetaboanalyst_annotatedDetected.svg',
    width=8, height=8, pointsize = 30)
venn(venn_annotated_detected, 
     ilabels = "counts", 
     ellipse = T, box = F, 
     zcolor = "style", 
     par = T)
dev.off()


annotated_detected %>% 
    filter(!pPTGE30) %>% 
    select(Clean_name, INCHIKEY, all_of(clones)) %>%
    write_delim(file = "Results/Annotated_notInE30.txt",
                delim="\t")
