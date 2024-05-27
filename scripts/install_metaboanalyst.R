#### MetaboAnalystR ####
# Step 1: Instalation of dependencies
# SSPA was removed in bioconductor version 3.14, install with url
install.packages("https://bioconductor.org/packages/3.12/bioc/bin/windows/contrib/4.0/SSPA_2.30.0.zip")
metanr_packages <- function(){
    metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "gtools",
                   "Rgraphviz", "preprocessCore", "genefilter", "SSPA", 
                   "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", 
                   "MSnbase", "multtest", "RBGL", "edgeR", "fgsea", "devtools", "crmn")
    
    list_installed <- installed.packages()
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install(new_pkgs)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}
metanr_packages()

# Step 2: Install MetaboAnalystR with documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, 
                         build_vignettes = TRUE, build_manual =T)
