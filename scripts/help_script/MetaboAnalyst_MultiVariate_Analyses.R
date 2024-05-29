# MetaboAnalyst multivariate analyses
# install.packages(c("ellipse", "pls"))
library(tidyverse)
library(MetaboAnalystR)
source("scripts/help_script/calculate.pls.vip.R")
theme_classic() %>% theme_set()

Analysis_modes = c("HILIC_Positive",
                   "RP_Positive",
                   "RP_Negative")
inFiles <- paste0("Results/", Analysis_modes, "_mSet.RData")
plot_basename <- paste0("plots/", Analysis_modes)
i = 1
for (i in 1:3) {
    rm("mSet")
    load(inFiles[i])
    mSet$dataSet$meta.info = 
        mSet$dataSet$meta.info %>% 
        filter(Class != "QC")
    mSet$dataSet$meta.info$Class <- 
        mSet$dataSet$meta.info$Class %>% 
        droplevels("QC")
    
    mSet$dataSet$cls = 
        mSet$dataSet$meta.info$Class 
    mSet$dataSet$norm = 
        mSet$dataSet$norm[1:nrow(mSet$dataSet$meta.info),]
    #### Perform PCA analysis ####
    mSet<-PCA.Anal(mSet)
    
    # Create PCA overview
    
    mSet<-PlotPCAPairSummary(
        mSet, pc.num = 3, 
        imgName = paste0(plot_basename[[i]], "_pca_pair_0_"),
        format = "png", dpi = 300, width=NA)
    
    # Create PCA scree plot
    mSet<-PlotPCAScree(
        mSet, 
        imgName = paste0(plot_basename[[i]], "_pca_scree_0_"),
        format = "png", dpi = 72, width=NA, 5)
    
    # Create a 2D PCA score plot
    mSet<-PlotPCA2DScore(
        mSet,
        imgName = paste0(plot_basename[[i]], "_pca_score2d_0_"), 
        format = "png", dpi=72, 
        width=NA, 1, 2, 0.95, 1, 0)
    

    # Create a PCA loadings Plots
    mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
    mSet$analSet$pca$imp.loads %>%
        as.data.frame() %>%
        mutate(Metabolite = rownames(.)) %>% 
        ggplot(aes(x = `Loadings 1`, `Loadings 2`)) +
            geom_point()
    
    ggsave(filename = paste0(plot_basename[[i]], "_pca_loading_0.png"),
           height = 5, width = 5, dpi = 300)

    #### Partial least square - discriminant analysis (PLS-DA) ####
    
    datmat <- as.matrix(mSet$dataSet$norm)
    cls <- model.matrix(~mSet$dataSet$cls - 1)
    pls.res <- pls::plsr(cls ~ datmat, method = "oscorespls", 
                         ncomp = 8, maxit = 10000)
    VIP_scores <- calculate.pls.vip(pls.res, 8)
    mSet$analSet$plsr <- pls.res
    mSet$analSet$plsr$reg <- F
    mSet$analSet$plsr$loading.type <- "all"
    mSet$analSet$plsr$vip.mat <- VIP_scores
    mSet$custom.cmpds <- c()
    
    plsda.cls <- caret::train(
        x=mSet$dataSet$norm, y=mSet$dataSet$cls, 
        "pls", maxit = 1000,
        trControl = caret::trainControl(method = "CV", 
                                        number = 5), 
        tuneLength = 8)
    
    plsda.reg <- pls::plsr(cls ~ datmat, 
                           method = "oscorespls", 
                           ncomp = 8, 
                           validation = "CV", maxit = 1000
                           )
    
    fit.info <- pls::R2(plsda.reg, estimate = "all")$val[, 1, ]
    accu <- plsda.cls$results[, 2]
    all.info <- rbind(accu, fit.info[, -1])
    rownames(all.info) <- c("Accuracy", "R2", "Q2")
    best.num <- min(which(all.info[2, ] == max(all.info[2, ])))
    
    coef.mat <- caret::varImp(plsda.cls, scale = T)$importance
    
    coef.mat <- cbind(coef.mean = rowMeans(coef.mat), 
                      coef.mat)
    
    inx.ord <- order(coef.mat[, 1], decreasing = T)
    coef.mat <- data.matrix(coef.mat[inx.ord, , drop = FALSE])
    mSet$analSet$plsda <- list(best.num = best.num, choice = "R2", 
                               coef.mat = coef.mat, fit.info = all.info)
    
    mSet<-PLSDA.CV(mSet, "5", compNum = 8, 
                   segments = 10, choice = "Q2")
    
    components = mSet$analSet$plsr$Xvar %>% 
        sort(decreasing = T) %>%
        names() %>% gsub(pattern="Comp ", replacement="") %>%
        as.numeric()
    
    mSet<-PlotPLSPairSummary(
        mSet, 
        imgName = paste0(plot_basename[[i]], "_pls_pair_0_"),
        "png", 72, width=NA, 5)
    
    mSet<-PlotPLS2DScore(
        mSet, 
        imgName = paste0(plot_basename[[i]], "_pls_score2d_0_"),
        "png", 72, width=NA, components[1],
        components[2],
        0.95,1,0)
    
    mSet<-PlotPLS3DScoreImg(
        mSet, 
        imgName = paste0(plot_basename[[i]], "_pls_score3d_0_"),
        "png", 72, width=NA, 
        components[1], components[2], components[3],
        40)

    
    mSet<-PlotPLS.Classification(
        mSet, 
        imgName = paste0(plot_basename[[i]], "_pls_cv_0_"),
        "png", 72, width=NA)
    
    mSet<-PlotPLS.Imp(mSet, 
                      imgName = paste0(plot_basename[[i]], "_pls_imp_0_"),
                      "png", 72, width=NA, "vip", "Comp. 1", 
                      10, FALSE)
    
    mSet<-PLSDA.Permut(mSet, 100, "accu")
}