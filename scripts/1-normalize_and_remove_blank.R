library(tidyverse)

blank_threshold = 3
in_files = list.files(path = "Inputs",
                      pattern = "^LCMSMS.*rawHeight.txt",
                      full.names = T)
out_files = sapply(in_files, \(x) {
    paste0("Inputs/CleanUp_", basename(x)) 
})

for (fileNumber in 1:3) {
    input = read_delim(in_files[fileNumber]) %>%
        as.data.frame()
    
    print(head(n=2, input))
    input[1,grepl("BLANK", input[1,])] = "BLANK"
    metadata = data.frame(Replicates = names(input)[-1],
                          Groups = input[1,-1] %>% as.character
    )
    
    measures = apply(input[-1,-1], 2, \(x) {
        as.numeric(x) + 1
    }) %>%
        as.data.frame
    rownames(measures) = input[-1,1]

    blank_value = apply(measures[, metadata$Groups == "BLANK"],
                        1, max)
    
    greater_than_blank = 
        apply(measures[, !metadata$Groups %in% c("BLANK", "QC")], 
              2, \(x) {
                  x > blank_value*blank_threshold
              }) %>% data.frame(row.names = rownames(measures))
    
    clones = metadata[!metadata$Groups %in% c("BLANK", "QC"),]$Groups %>% unique
    numberRepsPresent = 
        sapply(clones, \(x) {
            colsInterest = (metadata %>%
                                filter(Groups == x))$Replicates
            samples_GTB = greater_than_blank[, colsInterest] %>% rowSums() 
            samples_GTB > 1
        }) %>% as.data.frame
    
    measures_greater_than_blank = input %>%
        select(Sample,
               all_of(metadata$Replicates[metadata$Groups %in% clones])) %>%
        as.data.frame(row.names = input$Sample)
    
    measurements = sapply(names(measures_greater_than_blank)[-1],
           simplify = T, \(colName) {
            
            cloneName = metadata$Groups[metadata$Replicates == colName]
            sapply(1:nrow(numberRepsPresent), \(met) {
                ifelse(
                    numberRepsPresent[met, cloneName],
                    measures_greater_than_blank[met+1, colName],
                    0)
            })
        }) %>% as.data.frame(row.names = input$Sample[-1])
    measures = measures [!apply(measurements, 1, \(x) all(x == "0")),
                         metadata$Replicates[!metadata$Groups %in% c("BLANK")],]
    
    
    cbind(Sample = c("Label", rownames(measures)), 
          rbind((data.frame(metadata$Groups[metadata$Groups !="BLANK"], 
                     row.names = metadata$Replicates[metadata$Groups !="BLANK"]) %>% 
              unname %>% t), measures)
    ) %>% as_tibble %>%
        write_delim(file = out_files[fileNumber], delim = "\t", append = F)
}
