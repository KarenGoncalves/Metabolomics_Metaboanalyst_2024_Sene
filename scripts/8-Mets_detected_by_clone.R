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
nrow(venn_detected)

pdf('plots/Venn_filteredMetaboanalyst_detected.pdf',
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
nrow(venn_annotated_detected)

pdf('plots/Venn_filteredMetaboanalyst_annotatedDetected.pdf',
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

### Proportions ###
proportion_annotated = 
    Annotation %>%
    group_by(AnalysisMode) %>%
    summarize(`With annotation` = which(Clean_name != "Unknown") %>% length,
              Detected = n(),
              `Without annotation` = Detected - `With annotation`) %>%
    pivot_longer(cols = c(`With annotation`, `Without annotation`),
                 names_to = "Which",
                 values_to = "N") %>%
    mutate(Proportion = N/Detected,
           plotMode = gsub("RP_", "Reverse-phase\n", AnalysisMode) %>%
               gsub(pattern="_", replacement="\n"),
           Annotated = ifelse(Which == "With annotation", "Yes", "No")
           ) %>%
    select(plotMode, Annotated, Detected, N, Proportion)

write_delim("Results/Proportion_analytes_annotated.txt",
            x = proportion_annotated %>% select(!Proportion),
            quote = "none", delim="\t")

proportion_annotated %>%
    ggplot(aes(y = Proportion, x = "", fill = Annotated)) +
    geom_bar(stat = "identity") + 
    coord_polar(theta = 'y') +
    facet_grid(~plotMode, scales = "free_y",) +
    theme_void() +
    labs(fill = "") + 
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey70")) +
    theme(legend.position = "bottom",
          text=element_text(face = "bold", size = 10),
          strip.clip = "off")

ggsave("plots/Proportion_annotated_DAAs.pdf",
        height=2.5, width=5, dpi=1200)

