# MetaboAnalyst sigggenes pipeline
library(tidyverse)
library(MetaboAnalystR)
library(ggrepel)
library(RColorBrewer)
library(siggenes)
theme_classic() %>% theme_set()

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
inFiles <- paste0("Results/", Analysis_modes, "_mSet.RData")
plot_basename <- paste0("plots/", Analysis_modes, "/")
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("pPTGE30", "pPTGE30", "pPTGE30", 
                               "AC9.3", "AC9.3", "AC9.2"))
clones = c("AC9.1", "AC9.2", "AC9.3", "pPTGE30")
FDR_threshold = 0.01
FC_threshold = 1

met_abundance <- list()
differential_abundance_all <- list()
for (i in 1:3) {
    AM=Analysis_modes[i]
    rm("mSet")
    load(inFiles[i])
    
    met_abundance[[i]] <- mSet$dataSet$filt 
    differential_abundance <- list()
    
    for (j in 1:nrow(tableContrast)){
        print(tableContrast[j,])
        contrastName = paste0(tableContrast[j,1], "_v_", tableContrast[j,2])
        replicates = rownames(met_abundance[[i]])
        cls_subset = factor(rep(tableContrast[j,] %>% unlist, each = 3))
        
        if (tableContrast[j,2] == "pPTGE30")  {
            cls_subset <- relevel(cls_subset, "pPTGE30")
            cls_subset <- sort(cls_subset)
        }
        mSet_subset_norm = 
            met_abundance[[i]][c(grep(levels(cls_subset)[1], replicates),
                                 grep(levels(cls_subset)[2], replicates)), ] %>%
            t %>% as.data.frame()
        
        SAM = siggenes::sam(mSet_subset_norm, 
                            cl = cls_subset, 
                            method = "d.stat", B=200, 
                            gene.names = names(mSet_subset_norm), 
                            med = F, use.dm = T, 
                            var.equal = T, 
                            R.fold = 1
        )
        
        SAM@msg = c(contrastName, SAM@msg)
        
        FC <- sapply(rownames(mSet_subset_norm), \(x) {
            numerator_cols = grep(tableContrast[j,1], colnames(mSet_subset_norm))
            denominator_cols = grep(tableContrast[j,2], colnames(mSet_subset_norm))
            fold_change <- 
                (rowMeans(mSet_subset_norm[x,numerator_cols]) / 
                    rowMeans(mSet_subset_norm[x,denominator_cols])) %>% 
                log2
            ifelse(abs(fold_change) < FC_threshold,
                   0, fold_change)
        })
        
        differential_abundance[[contrastName]] <- 
            data.frame(AnalysisMode = AM,
                       Contrast = contrastName,
                       plotContrast = paste0(tableContrast[j,1], " vs ", tableContrast[j,2]),
                       Metabolite = rownames(mSet_subset_norm),
                       Fold = 1/SAM@fold,
                       pValue = SAM@p.value,
                       padj = SAM@q.value)
        differential_abundance[[contrastName]]$FoldChange = 
            sapply(1:nrow(differential_abundance[[contrastName]]), \(x) {           
                with(differential_abundance[[contrastName]],
                     ifelse(pValue[x] < FDR_threshold, FC[x], 0)
                )
            })
    }
    differential_abundance_all[[AM]] <- 
        differential_abundance %>%
        list_rbind()
}

for (i in list.files(".", patter=".(qs|csv)")) {
    # remove these temporary files
    file.remove(i)
}

differential_abundance_all %>% list_rbind() %>%
    write_delim(file = "Results/Siggenes_DAAs.txt",
                quote="none", append = F, delim="\t")
