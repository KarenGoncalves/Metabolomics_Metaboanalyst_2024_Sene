#'Comparison of means, unpaired
#'@description performs homoscedascity test followed by either ANOVA+TukeyHSD or Kruskal-Wallis+Dunn's test
#'@usage ANOVA_mets(data, clones, mets)
#'@param data mSetObj normalized object (eg. mSet$dataSet$norm[1:nrow(clones), ] %>% as.data.frame)
#'@param mets names of metabolites (column names of data)
#'@param clones table of indicating to which treatment each row belongs to (nrow(clones) == nrow(data))
#'@author Karen Goncalves\email{cris.kgs@gmail.com}
#'UQTR, Canada
#'@import tidyverse
#'@import car
#'@import dunn.test
#'@export

ANOVA_mets <- function(data, clones, mets) {
    require(dunn.test)
    require(car)
    require(tidyverse)
    
    result = sapply(mets, simplify = F, \(x) {
        test <- list()
        nonpar = leveneTest(data[[x]] ~ 
                                clones[[1]])$`Pr(>F)`[1]
        test[["homoscedascity"]] <- nonpar
        if (nonpar < .05) {
            dif_test <- 
                kruskal.test(data[[x]] ~ 
                                 clones[[1]])
            
            if (dif_test$p.value < .05) {
                posthoc_test <-  
                    dunn.test(x = data[[x]],
                              g = clones[[1]]) %>% as.data.frame 
                test[["posthoc"]] <- posthoc_test %>%
                    separate(comparisons, into = c("Clone1", "Clone2"), sep = " - ") %>%
                    rename("pvalue_posthoc" = P) %>%
                    mutate(Metabolite = x,
                           Method = dif_test$method,
                           statistic = dif_test$statistic %>% unname,
                           Posthoc = "Dunn's test",
                           ANOVA_p.value = dif_test$p.value,
                           ANOVA_padj = dif_test$p.value %>% 
                               p.adjust(n = length(mets))
                    ) %>%
                    select(Metabolite, Method, Posthoc, statistic, ANOVA_p.value, 
                           ANOVA_padj, Clone1, Clone2, pvalue_posthoc)
                
            } else {
                test[["posthoc"]] <- 
                    data.frame(Clone1 = tableContrast$Numerator,
                               Clone2 = tableContrast$Denominator,
                               Metabolite = x,
                               Method = dif_test$method,
                               statistic = dif_test$statistic %>% unname,
                               Posthoc = NA,
                               ANOVA_p.value = dif_test$p.value,
                               ANOVA_padj = dif_test$p.value %>% 
                                   p.adjust(n = length(mets)),
                               pvalue_posthoc = 1
                    ) %>%
                    select(Metabolite, Method, Posthoc, statistic, ANOVA_p.value, 
                           ANOVA_padj, Clone1, Clone2, pvalue_posthoc)
            }
        } else {
            dif_test = 
                aov(data[[x]] ~ clones[[1]])
            pval <- (dif_test%>% 
                         summary)[[1]][[5]][[1]] 
            if (pval < .05) {
                test[["posthoc"]] <- 
                    pairwise.t.test(x = data[[x]], 
                                    g = clones[[1]], 
                                    paired = F)[["p.value"]] %>% 
                    as.data.frame() %>% mutate(Clone1 = row.names(.)) %>%
                    pivot_longer(cols = !Clone1,
                                 names_to =  "Clone2",
                                 values_to = "pvalue_posthoc", values_drop_na = T) %>%
                    mutate(Metabolite = x,
                           Method = "ANOVA",
                           statistic = (dif_test %>% summary())[[1]][[4]][[1]],
                           Posthoc = "Tukey's test",
                           ANOVA_p.value = pval,
                           ANOVA_padj = p.adjust(p = pval, n = length(mets))
                    ) %>%
                    select(Metabolite, Method, Posthoc, statistic, ANOVA_p.value, 
                           ANOVA_padj, Clone1, Clone2, pvalue_posthoc)
            } else {
                
                test[["posthoc"]] <- 
                    data.frame(Clone1 = tableContrast$Numerator,
                               Clone2 = tableContrast$Denominator,
                               Metabolite = x,
                               Method = "ANOVA",
                               statistic = (dif_test %>% summary())[[1]][[4]][[1]],
                               Posthoc = NA,
                               ANOVA_p.value = pval,
                               ANOVA_padj = p.adjust(p = pval, n = length(mets)),
                               pvalue_posthoc = 1
                    ) %>%
                    select(Metabolite, Method, Posthoc, statistic, ANOVA_p.value, 
                           ANOVA_padj, Clone1, Clone2, pvalue_posthoc)
            }
        }
        test
    })
    return(result)
}
