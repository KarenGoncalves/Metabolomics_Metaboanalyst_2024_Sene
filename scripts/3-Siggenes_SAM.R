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
FDR_threshold = 0.05
adjusted_FDR_threshold = 0.00833333233333
FC_threshold = 2

met_abundance <- list()
differential_abundance_all <- list()
for (i in 1:3) {
    AM=Analysis_modes[i]
    pdf(file = paste0(plot_basename[i], 
                      "Differential_analysis.pdf"),
        width = 6, height = 6, onefile = T)
    rm("mSet")
    load(inFiles[i])
    
    met_abundance[[i]] <- mSet$dataSet$norm 
    differential_abundance <- list()
    
    for (j in 1:nrow(tableContrast)){
        print(tableContrast[j,])
        contrastName = paste0(tableContrast[j,1], "_v_", tableContrast[j,2])
        replicates = rownames(mSet$dataSet$norm)
        cls_subset = factor(rep(tableContrast[j,] %>% unlist, each = 3))
        
        if (tableContrast[j,2] == "pPTGE30")  {
            cls_subset <- relevel(cls_subset, "pPTGE30")
            cls_subset <- sort(cls_subset)
        }
        mSet_subset_norm = 
            met_abundance[[i]][c(grep(levels(cls_subset)[1], replicates),
                                 grep(levels(cls_subset)[2], replicates)), ]
        SAM = siggenes::sam(t(mSet_subset_norm), 
                            cl = cls_subset, 
                            method = "d.stat", B=100, 
                            gene.names = names(mSet_subset_norm), 
                            med = F, use.dm = T, 
                            var.equal = T,
                            control = samControl(n.delta = 120),
                            R.fold = FC_threshold
                            
        )
        SAM@msg = c(contrastName, SAM@msg)
        
        deltaValue = findDelta(SAM, fdr = adjusted_FDR_threshold)
        if (nrow(deltaValue) %>% is.null) {
            deltaValue = deltaValue["Delta"]
        } else {
            deltaValue = deltaValue[order(deltaValue[,"FDR"]), "Delta"]
        }
        sam.plot2(SAM, deltaValue, sig.col = c("blue", "red"), 
                  main = paste0(tableContrast[j,1], " vs ", tableContrast[j,2])
        )
        
        differential_abundance[[contrastName]] <- 
            data.frame(AnalysisMode = AM,
                       Contrast = contrastName,
                       plotContrast = paste0(tableContrast[j,1], " vs ", tableContrast[j,2]),
                       Metabolite = names(SAM@d),
                       Ratio_Change = 1/SAM@fold,
                       FoldChange = log2(1/SAM@fold),
                       pValue = SAM@p.value,
                       padj = SAM@q.value) 
                         
    }
    differential_abundance_all[[AM]] <- 
        differential_abundance %>%
        list_rbind()
dev.off()
}

for (i in list.files(".", patter=".(qs|csv)")) {
    file.remove(i)
}

save(differential_abundance_all, file = "Results/Siggenes_DAAs.RData")
