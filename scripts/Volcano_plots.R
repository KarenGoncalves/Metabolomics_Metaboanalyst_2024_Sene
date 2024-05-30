# Plots differential regulation
library(tidyverse)
library(MetaboAnalystR)

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
tableContrast <- 
    data.frame(Numerator = c("AC9.1", "AC9.2", "AC9.3",
                             "AC9.1", "AC9.2", "AC9.1"),
               Denominator = c("E30", "E30", "E30", 
                               "AC9.3", "AC9.3", "AC9.2"))

DAM_tables <- sapply(Analysis_modes, simplify = F, \(x) {
    read_delim(paste0("Results/", x, "_DAM_table.txt")) %>%
        mutate(AnalysisMode = x)
}) %>% list_rbind


FC_res_long <- sapply(Analysis_modes, simplify = F, \(x) {
    read_delim(paste0("Results/", x, "_FC_result.txt")) %>%
        mutate(AnalysisMode = x)
}) %>% list_rbind

inverse_contrast <- function(tableContrast, DAM_table) {
    new_DAM_table = DAM_tables
    for (i in 1:nrow(tableContrast)) {
        reverse_rows = which(DAM_table$Clone1 == tableContrast$Denominator[i] &
                                 DAM_table$Clone2 == tableContrast$Numerator[i])
        new_DAM_table[reverse_rows, c("Clone1", "Clone2")] =
            tableContrast[i,]
    }
    return(new_DAM_table)
}

FC_ANOVA <- FC_res_long %>% 
    full_join(inverse_contrast(tableContrast, DAM_tables),
              by = c("Clone1", "Clone2", "Metabolite", "AnalysisMode")) %>%
    unique %>%
    mutate(Deregulation = 
               case_when(log2FC > 2 & pvalue_posthoc < .05 ~ "Up-regulated",
                         log2FC < -2 & pvalue_posthoc < .05 ~ "Down-regulated",
                         .default = "No"))
write_delim(FC_ANOVA, 
            "Results/FC_ANOVA_allmodes.txt",
            delim = "\t")
FC_limit = max(abs(FC_ANOVA$log2FC))

FC_ANOVA$plotContrast <- 
    FC_ANOVA$plotContrast %>% 
    factor(levels = 
               apply(tableContrast, 1, \(x) paste(x[1], x[2], 
                                                  sep = " vs "))
    )

Number_deregulated <- FC_ANOVA %>%
    group_by(plotContrast, Deregulation, AnalysisMode) %>%
    summarize(Total = length(Metabolite)) %>%
    filter(Deregulation != "No") %>%
    mutate(x = case_when(Deregulation == "Up-regulated" ~ 5,
                         .default = -5),
           y = 10)

for (i in Analysis_modes) {
    FC_ANOVA %>%
        filter(pvalue_posthoc != 1 &
                   AnalysisMode == i) %>%
        ggplot(aes(x = log2FC, y = -log10(pvalue_posthoc),
                   color = Deregulation, alpha = Deregulation)) +
        geom_point(stroke=.3) +
        geom_hline(yintercept = -log10(.05), 
                   linetype = 2, color = "grey") +
        geom_vline(xintercept = c(-2, 2), 
                   linetype = 2) +
        xlim(c(-FC_limit, FC_limit)) +
        facet_wrap(~plotContrast, 
                   nrow = 2, ncol = 3, as.table = T) +
        geom_text(data = Number_deregulated %>%
                      filter(AnalysisMode == i),
                  aes(x, y = y, label = Total),
                  show.legend = F
        ) +
        scale_color_manual(values = c("Down-regulated" = "blue",
                                      "Up-regulated"="red",
                                      "No" = "grey50"),
                           name="") +
        scale_alpha_manual(values = c("Up-regulated" = 1,
                                      "Down-regulated" = 1,
                                      "No" = .5),
                           name="") +
        labs(x = bquote("log"[2]~"Fold Change"),
             y = bquote("-log"[10]~"pValue")) +
        theme_classic() +
        theme(legend.position = "bottom") 
    ggsave(paste0("plots/", i, "/Volcano.png"),
           height=5, width=8)
}

#### Volcano combined

Number_deregulated <- FC_ANOVA %>%
    group_by(plotContrast, Deregulation) %>%
    summarize(Total = length(Metabolite)) %>%
    filter(Deregulation != "No") %>%
    mutate(x = case_when(Deregulation == "Up-regulated" ~ 5,
                         .default = -5),
           y = 10)

FC_ANOVA %>%
    filter(pvalue_posthoc != 1 ) %>%
    ggplot(aes(x = log2FC, y = -log10(pvalue_posthoc),
               color = Deregulation, alpha = Deregulation)) +
    geom_point(stroke=.3) +
    geom_hline(yintercept = -log10(.05), 
               linetype = 2, color = "grey") +
    geom_vline(xintercept = c(-2, 2), 
               linetype = 2) +
    xlim(c(-FC_limit, FC_limit)) +
    facet_wrap(~plotContrast, 
               nrow = 2, ncol = 3, as.table = T) +
    geom_text(data = Number_deregulated,
              aes(x, y = y, label = Total),
              show.legend = F
    ) +
    scale_color_manual(values = c("Down-regulated" = "blue",
                                  "Up-regulated"="red",
                                  "No" = "grey50"),
                       name="") +
    scale_alpha_manual(values = c("Up-regulated" = 1,
                                  "Down-regulated" = 1,
                                  "No" = .5),
                       name="") +
    labs(x = bquote("log"[2]~"Fold Change"),
         y = bquote("-log"[10]~"pValue")) +
    theme_classic() +
    theme(legend.position = "bottom") 
ggsave(paste0("plots/Volcano_allModes.png"),
       height=5, width=8)



