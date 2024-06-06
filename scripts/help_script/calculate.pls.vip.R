#'Get metabolite importance from pls
#'@description returns the VIP scores of each metabolite
#'@usage calculate.pls.vip(plsr.result, compNum) 
#'@param data result from pls::plsr function
#'@param compNum number of components used in plsr function
#'@export
calculate.pls.vip <- function(plsr.result, compNum) {
    pls <- plsr.result
    b <- c(pls$Yloadings)[1:compNum]
    T <- pls$scores[, 1:compNum, drop = FALSE]
    SS <- b^2 * colSums(T^2)
    W <- pls$loading.weights[, 1:compNum, drop = FALSE]
    Wnorm2 <- colSums(W^2)
    SSW <- sweep(W^2, 2, SS/Wnorm2, "*")
    vips <- sqrt(nrow(SSW) * apply(SSW, 1, cumsum)/cumsum(SS))
    if (compNum > 1) {
        vip.mat <- as.matrix(t(vips))
        ord.inx <- order(-abs(vip.mat[, 1]), -abs(vip.mat[, 
                                                          2]))
    }
    else {
        vip.mat <- as.matrix(vips)
        ord.inx <- order(-abs(vip.mat[, 1]))
    }
    vip.mat <- vip.mat[ord.inx, ]
    colnames(vip.mat) <- paste("Comp.", 1:ncol(vip.mat))
    return(vip.mat)
}

#' Analyse and plot PLS-DA
#'@description returns the VIP scores of each metabolite
#'@usage runPLS.Analysis(mSetObj, compNum=10, iterations=1000, method = "oscorespls", basenamePlots) 
#'@param data result from pls::plsr function
#'@param compNum number of components used in plsr function
#'@import pls
#'@import caret
#'@import tidyverse
#'@import MetaboAnalyst
#'@export
runPLS.Analysis <- function(mSetObj, compNum, iterations, 
                            method = "oscorespls", basenamePlots) {
    require(pls)
    require(caret)
    require(tidyverse)
    require(MetaboAnalystR)
    datmat <- as.matrix(mSetObj$dataSet$norm)
    cls <- model.matrix(~mSetObj$dataSet$cls - 1)
    pls.res <- pls::plsr(cls ~ datmat, method = method, 
                         ncomp = compNum, maxit = iterations)
    VIP_scores <- calculate.pls.vip(pls.res, compNum = compNum)
    mSetObj$analSet$plsr <- pls.res
    mSetObj$analSet$plsr$reg <- F
    mSetObj$analSet$plsr$loading.type <- "all"
    mSetObj$analSet$plsr$vip.mat <- VIP_scores
    mSetObj$custom.cmpds <- c()
    
    plsda.cls <- caret::train(
        x=mSetObj$dataSet$norm, y=mSetObj$dataSet$cls, 
        "pls", maxit = iterations,
        trControl = caret::trainControl(method = "CV", 
                                        number = 5), 
        tuneLength = compNum)
    
    plsda.reg <- pls::plsr(cls ~ datmat, 
                           method = method, 
                           ncomp = compNum,
                           validation = "CV", 
                           maxit = iterations
    )
    
    fit.info <- pls::R2(plsda.reg, estimate = "all", )$val[, 1, ]
    
    accu <- plsda.cls$results[, 2]
    all.info <- rbind(accu, fit.info[, -1])
    rownames(all.info) <- c("Accuracy", "R2", "Q2")
    best.num <- min(which(all.info[2, ] == max(all.info[2, ])))
    
    coef.mat <- caret::varImp(plsda.cls, scale = T)$importance
    
    coef.mat <- cbind(coef.mean = rowMeans(coef.mat), 
                      coef.mat)
    
    inx.ord <- order(coef.mat[, 1], decreasing = T)
    coef.mat <- data.matrix(coef.mat[inx.ord, , drop = FALSE])
    mSetObj$analSet$plsda <- list(best.num = best.num, choice = "R2", 
                                  coef.mat = coef.mat, fit.info = all.info)
    
    mSetObj<-PlotPLSPairSummary(
        mSetObj,
        imgName = paste0(basenamePlots, "PLS_pair_0_"),
        "png", 72, width=NA, 5)
     
    mSetObj<-PlotPLS2DScore(
        mSetObj,
        imgName = paste0(basenamePlots, "PLS_score2d_0_"),
        "png", 72, width=NA, inx1 = 1, inx2 = 2,
        reg = 0.95,show = 1,grey.scale = 0)
    
    mSetObj<-PlotPLS.Classification(
        mSetObj,
        imgName = paste0(basenamePlots, "PLS_cv_0_"),
        "png", 72, width=NA)
     
    mSetObj<-PlotPLS.Imp(
        mSetObj,
        imgName = paste0(basenamePlots, "PLS_imp_0_"),
        "png", 72, width=NA, type = "vip", feat.nm = "Comp. 1",
        feat.num=10, FALSE)
     
    return(mSetObj)
}

plot_RegCoef <- \(mSetObj, comp){
    pls.res = mSetObj$analSet$plsr
    coef_plot = pls.res$coefficients[, , comp] %>%
        as.data.frame() 
    names(coef_plot) = gsub("mSetObj\\$dataSet\\$cls", "", names(coef_plot))
    coef_plot = coef_plot %>% mutate(Metabolites =  rownames(.)) %>%
        pivot_longer(cols = !Metabolites,
                     names_to = "Clones",
                     values_to = "Regression_coef") %>%
        ggplot(aes(x = Metabolites, y = Regression_coef)) +
        geom_col() +
        geom_hline(yintercept=0) +
        labs(title = paste("Component", comp),
             y = "Regression coefficient") +
        facet_wrap(~Clones, ncol = 2) +
        theme_classic() +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
    return(coef_plot)
}
