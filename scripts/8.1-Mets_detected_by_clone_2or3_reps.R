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

featuresKept_metaboanalyst <- sapply(mSetData, \(x) {
    load(x) 
    mSet[["dataSet"]][["filt"]] %>% colnames
})%>% unlist %>% unname %>% c
 
input = sapply(1:length(Analysis_modes), simplify = F, \(i) {
    read_delim(inFiles[i])[-1,] 
}) %>% list_rbind() %>%
    select(all_of(c("Sample", 
                    paste0(rep(clones, each=3), 
                           "_", c("I", "II", "III")
                    ))
    )) %>% 
    rename(Metabolite = Sample) %>%
    filter(Metabolite %in% featuresKept_metaboanalyst)

analytes_detected <- 
    sapply(clones, simplify = T, \(x) {
        data_clone <- 
            input %>%
            select(
                all_of(paste0(x, "_", c("I", "II", "III")))
            )
        apply(data_clone, 2, \(y) as.numeric(y) %>% as.logical) %>%
            rowSums
    }) %>% as_tibble() 

#### Detected in at least 2 reps ####
analytes_detected_2orMore <- 
    apply(analytes_detected, 2, \(x) {
        x >= 2
        })  %>%
    as.data.frame(row.names = input$Metabolite) %>%
    filter(rowSums(.) != 0) %>%
    mutate(Metabolite = row.names(.)) %>%
    separate(col = Metabolite, into = c("RT", "Mz"), 
             sep = "/", convert = T)


Annotation_2orMore <- 
    read_delim("Results/Annotation_clean_names.txt") %>%
    right_join(analytes_detected_2orMore, 
               by = join_by("Average Rt(min)" == "RT",
                            "Average Mz" == "Mz")) 

annotated_detected_2orMore <- Annotation_2orMore %>%
    filter(Clean_name != "Unknown" &
               !grepl("w/o MS2", Clean_name))



#### Detected in all reps ####

analytes_detected_all3 <- 
    apply(analytes_detected, 2, \(x) {
        x == 3
    }) %>%
    as.data.frame(row.names = input$Metabolite) %>%
    filter(rowSums(.) != 0) %>%
    mutate(Metabolite = row.names(.)) %>%
    separate(col = Metabolite, into = c("RT", "Mz"), 
             sep = "/", convert = T)

Annotation_all3 <- 
    read_delim("Results/Annotation_clean_names.txt") %>%
    right_join(analytes_detected_all3, 
               by = join_by("Average Rt(min)" == "RT",
                            "Average Mz" == "Mz")) 

annotated_detected_all3 <- 
    Annotation_all3 %>%
    filter(Clean_name != "Unknown" &
               !grepl("w/o MS2", Clean_name))

#### Venn plots ####
##### 2 or more #####
par(cex=1, font=1); 

pdf('plots/Venn_filteredMetaboanalyst_2orMore.pdf',
    width=8, height=8, pointsize = 30)
venn(analytes_detected_2orMore[, clones], main="Detected",
     ilabels = "counts", 
     ellipse = T, box = F, 
     zcolor = "style", 
     par = T)
venn(annotated_detected_2orMore[, clones], main="Annotated", 
     ilabels = "counts", 
     ellipse = T, box = F, 
     zcolor = "style", 
     par = T)

dev.off()

##### All three reps #####
pdf('plots/Venn_filteredMetaboanalyst_all3.pdf',
    width=8, height=8, pointsize = 30)
venn(analytes_detected_all3[, clones], main="Detected",
     ilabels = "counts", 
     ellipse = T, box = F, 
     zcolor = "style", 
     par = T)
venn(annotated_detected_all3[, clones], main="Annotated", 
     ilabels = "counts", 
     ellipse = T, box = F, 
     zcolor = "style", 
     par = T)

dev.off()

#### table ####
annotated_detected_2orMore %>% 
    filter(!pPTGE30) %>% 
    select(Clean_name, INCHIKEY, all_of(clones)) 

annotated_detected_all3 %>% 
    filter(!pPTGE30) %>% 
    select(Clean_name, INCHIKEY, all_of(clones)) 
