# MetaboAnalyst Analysis
install.packages("agricolae")
library(tidyverse)
library(MetaboAnalystR)
library(car)
library(dunn.test)
library(agricolae)


i=1
in_files = list.files(path = "Inputs/",
                      pattern = "CleanUp_",
                      full.names = T)
AnalysisMode <- gsub("CleanUp_LCMSMS_(.+)_rawHeight.txt",
                     "\\1", basename(in_files[1]))
plots_basename <- paste0("plots/", AnalysisMode)

mSet<-InitDataObjects("conc", "stat", FALSE);
mSet<-Read.TextData(mSet, in_files[1], "colu", "disc");
mSet<-SanityCheckData(mSet);
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet);
mSet<-FilterVariable(mSet, qc.filter = "T", rsd = 25, 
                     var.filter = "none", var.cutoff = -1, 
                     int.filter = "mean", int.cutoff = 0)
mSet<-PreparePrenormData(mSet);
mSet<-Normalization(mSet, "GroupPQN", transNorm = "LogNorm", 
                    scaleNorm = "ParetoNorm", ref = "QC", 
                    ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, imgName = paste0(plots_basename,
                                             "_rowNormalization_"), 
                      format ="png", dpi=300, width=NA);
mSet<-PlotSampleNormSummary(mSet, imgName = paste0(plots_basename,
                                                   "_colNormalization_"), 
                            format ="png", dpi=300, width=NA);

# Perform fold-change analysis on uploaded data, unpaired
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("E30", "E30", "E30", 
                               "AC9.3", "AC9.3", "AC9.2"))
FC_res<- MetaboAnalystR::GetFC(mSetObj = mSet, paired = F, 
                             tableContrast = tableContrast)

# Plot fold-change analysis
mSet$analSet$fc <- FC_res
FC_res_table <- data.frame(IDs = FC_res$AC9.1_v_E30$fc.log %>% names)
for (i in names(FC_res)) {
    FC_res_table[[i]] = FC_res[[i]]$fc.log
}

FC_res_table %>% pivot_longer(cols = !IDs,
                              values_to = "log2FC",
                              names_to = "Contrasts") %>%
    mutate(plotContrast = gsub("_v_", " vs ", Contrasts),
           Regulation = case_when(log2FC > 2 ~ "Up-regulated",
                                  log2FC < -2 ~ "Down-regulated",
                                  .default = "No")) %>%
    ggplot(aes(y = log2FC, x = IDs, 
               color = Regulation, alpha = Regulation)) +
    geom_point() + 
    scale_color_manual(values = c("Up-regulated" = "red",
                                  "Down-regulated" = "blue",
                                  "No" = "grey50"),
                       name = "") +
    scale_alpha_manual(values = c("Up-regulated" = 1,
                                  "Down-regulated" = 1,
                                  "No" = 0.3), 
                       name = "") +
    facet_wrap(~plotContrast, nrow = 3, ncol = 2, as.table = F) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    labs(x = "Analytes", y = bquote("log"[2]~"FC"))
ggsave(paste0(plots_basename, "_FCAnalysis.png"), dpi = 300,
       height = 7, width = 6)

# To view fold-change 
mets <- names(mSet$dataSet$norm)
clones <- mSet$dataSet$meta.info %>%
    filter(Class != "QC")
data <- mSet$dataSet$norm[1:nrow(clones), ] %>% as.data.frame()
anova_result <- 
    sapply(mets, simplify = F, \(x) {
        test <- list()
        nonpar = leveneTest(data[[x]] ~ 
                                clones[[1]])$`Pr(>F)`[1]
        test[["homoscedascity"]] <- nonpar
        if (nonpar < .05) {
            test[["anova"]] = 
                kruskal.test(data[[x]] ~ 
                                 clones[[1]])
            if (test[["anova"]]$p.value < .05) {
                test[["posthoc"]] <- 
                    dunn.test(x = data[[x]],
                              g = clones[[1]]) %>% as.data.frame
                test[["posthoc"]]$P.adjusted <- 
                    p.adjust( test[["posthoc"]]$P, n = length(mets))
            } else {
                test[["posthoc"]] <- NULL
            }
        } else {
            test[["anova"]] = 
                aov(data[[x]] ~ clones[[1]])
            pval <- (test[["anova"]]%>% 
                summary)[[1]][[5]][[1]]
            if (pval < .05) {
                test[["posthoc"]] <- 
                    pairwise.t.test(x = data[[x]], 
                                    g = clones[[1]], 
                                    paired = F)[["p.value"]] %>% 
                    as.data.frame() %>% mutate(Clone = row.names(.)) %>%
                    pivot_longer(cols = !Clone,
                                 names_to =  "Clone2",
                                 values_to = "p.value", values_drop_na = T) %>%
                    mutate(padj = p.adjust(p.value, 
                                           n = length(mets)))
                    
            } else {
                test[["posthoc"]] <- NULL
            }
        }
        test
    })

