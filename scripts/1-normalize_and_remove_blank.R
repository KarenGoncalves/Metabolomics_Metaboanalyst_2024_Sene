library(tidyverse)

blank_threshold = 3
in_files = list.files(path = "Inputs",
                      pattern = "^LCMSMS.*rawHeight.txt",
                      full.names = T)
out_files = sapply(in_files, \(x) {
    paste0("Inputs/CleanUp_", basename(x)) 
})
new_emptyVector = "pPTGE30"
min_reps_pass_threshold=3

for (fileNumber in 1:3) {
    input = read_delim(in_files[fileNumber]) %>%
        as.data.frame()
    names(input) <- gsub("E30", new_emptyVector, names(input))
    print(head(n=2, input))
    input[1,grepl("BLANK", input[1,])] = "BLANK"
    metadata = data.frame(Replicates = names(input)[-1],
                          Groups = input[1,-1] %>% as.character
    )
    
    measures = apply(input[-1,-1], 2, \(x) {
        as.numeric(x)
    }) %>%
        as.data.frame
    rownames(measures) = input[-1,1]

    blank_value = apply(measures[, metadata$Groups == "BLANK"],
                        1, max)
    
    clones = metadata[!metadata$Groups %in% c("BLANK", "QC"),]$Groups %>% unique
    greater_than_blank = 
        apply(measures[, !metadata$Groups  == "BLANK"], 
              2, \(x) {
                  result = (x - blank_value*blank_threshold)
                  ifelse(result < 0, 0,
                         result)
              }) %>% data.frame(row.names = rownames(measures))
    
    present_twoMore_reps = 
        sapply(c(clones, "QC"), \(x) {
            colsInterest = (metadata %>%
                                filter(Groups == x))$Replicates
            if (min_reps_pass_threshold < length(colsInterest)) {
                samples_GTB = greater_than_blank[, colsInterest] %>% 
                    apply(MARGIN=1, \(x) as.logical(x) %>% sum)
                samples_GTB >= min_reps_pass_threshold    
            } else {
                greater_than_blank[, colsInterest] %>% 
                    apply(MARGIN=1, \(x) as.logical(x) %>% all)
            }
        }) %>% as.data.frame 
    
    measures_greater_than_blank = greater_than_blank %>%
        mutate(Sample = input$Sample[-1]) %>%
        select(Sample,
               all_of(metadata$Replicates[metadata$Groups != "BLANK"])) %>%
        as.data.frame(row.names = input$Sample[-1])
    
    measurements = sapply(names(measures_greater_than_blank)[-1],
           simplify = T, \(colName) {
            
            cloneName = metadata$Groups[metadata$Replicates == colName]
            sapply(1:nrow(present_twoMore_reps), \(met) {
                ifelse(
                    present_twoMore_reps[met, cloneName],
                    measures_greater_than_blank[met, colName],
                    0)
            })
        }) %>% as.data.frame(row.names = input$Sample[-1])
    
    measures = 
        measurements[
            !apply(measurements[,metadata$Replicates[metadata$Groups %in% clones]],
                   1, \(x) all(x == "0")),
            metadata$Replicates[!metadata$Groups %in% c("BLANK")],]
    
    
    cbind(Sample = c("Label", rownames(measures)), 
          rbind((data.frame(metadata$Groups[metadata$Groups !="BLANK"], 
                     row.names = metadata$Replicates[metadata$Groups !="BLANK"]) %>% 
              unname %>% t), measures)
    ) %>% as_tibble %>%
        write_delim(file = out_files[fileNumber], delim = "\t", append = F)
}
