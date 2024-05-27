# Adapting MetaboAnalystR functions
library(ggplot2)
library(scales)
mSetObj$imgSet$fc <- list()
PlotFC_multiple <- 
function (mSetObj = NA, ContrastName, fc.thresh,
          imgName, format = "png", dpi = 72, width = NA) 
{
    imgName = paste(imgName, "dpi", dpi, ".", format, sep = "")
    if (is.na(width)) {
        w <- 8
    } else if (width == 0) {
        w <- 7
    } else {
        w <- width
    }
    h <- w * 6/8
    mSetObj$imgSet$fc[[ContrastName]] <- imgName
    fc = mSetObj$analSet$fc[[ContrastName]]
    fc_data <- data.frame(x = 1:length(fc$fc.log), y = fc$fc.log, 
                          Significance = ifelse(fc$inx.imp, "Significant", "Unsignificant"), 
                          label = names(fc$inx.imp))
    topVal <- max(abs(fc$fc.log))
    ylim <- c(-topVal, topVal)
    fc_data$ColorValue <- ifelse(fc_data$Significance == "Significant", 
                                 fc_data$y, NA)
    fc_data$ColorValue <- as.numeric(fc_data$ColorValue)
    sig_fc_range <- range(fc_data$ColorValue, na.rm = TRUE)
    fc_data$tooltip <- paste0("Label: ", fc_data$label, "\nLog2FC: ", 
                              round(fc_data$y, 2))
    nUp = filter(fc_data, fc.log > fc.thresh & Significance == "Significant") %>% nrow
    nDown = filter(fc_data, fc.log < -fc.thresh & Significance == "Significant") %>% nrow
    p <- ggplot(fc_data, aes(x = x, y = y, label = label, text = tooltip)) + 
        geom_point(aes(size = abs(y), color = ColorValue, fill = ColorValue), 
                   shape = 21, stroke = 0.5) + 
        scale_color_gradient2(low = "blue", 
                              mid = "grey", high = "red", midpoint = 0, limits = sig_fc_range, 
                              space = "Lab", na.value = "darkgrey", guide = "colourbar", 
                              name = "Log2FC") + 
        scale_fill_gradient2(low = "blue", 
                             mid = "grey", high = "red", midpoint = 0, limits = sig_fc_range, 
                             space = "Lab", na.value = "darkgrey", name = "Log2FC") + 
        scale_size(range = c(1.5, 4), guide = "none") + 
        labs(x = "Identifier", 
             y = "Log2 Fold Change",
             caption = paste0("Metabolites deregulated in ", 
                              gsub("_v_", " vs ", ContrastName),
                              ": Up=", nUp, "; Down=", nDown)) + 
        theme_bw() + theme(legend.position = "right") + 
        geom_hline(yintercept = 0, linewidth = 1)
        Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, 
                     width = w, height = h, type = format, bg = "white")
        print(p)
        dev.off()
        return(mSetObj)
}


# ANOVA

# aof <- function (x, cls) {
#     require(car)
#     if (leveneTest(x ~ cls)$`Pr(>F)`[1] > .05) {
#         aov(x ~ cls)
#     } else {
#         kruskal.test(x ~ cls)
#     }
# }
# 
# get.ftest.res <- function (data, cls) {
#     cls = cls[cls != "QC"]
#     print("Performing regular ANOVA F-tests ....")
#     aov.res <- my.res <-  NULL
#         aov.res <- apply(data, 2, aof, cls)
#         anova.res <- lapply(aov.res, anova)
#         my.res <- unlist(lapply(anova.res, function(x) {
#             c(x["F value"][1, ], x["Pr(>F)"][1, ])
#         }))
#     }
#     else {
#         require(dunn.test)
#         anova.res <- apply(data, 2, kwtest, cls)
#         my.res <- unlist(lapply(anova.res, function(x) {
#             c(x$statistic, x$p.value)
#         }))
#     }
#     my.res <- list(aov.res = aov.res, 
#                    f.res = data.frame(matrix(my.res, 
#                                              nrow = ncol(data), byrow = T), 
#                                       stringsAsFactors = FALSE))
#     return(my.res)
# }
GetFtestRes <- function(mSetObj=NA, nonpar=F){
    
    if(!exists("mem.aov")){
        require("memoise");
        mem.aov <<- memoise(get.ftest.res);
    }
    
    # mSetObj <- .get.mSet(mSetObj);  
    data <- as.matrix(mSetObj$dataSet$norm);
    cls <- mSetObj$dataSet$cls;
    print("using cache .......");
    return(mem.aov(data, cls, nonpar));
}

# ANOVA.Anal <- function (mSetObj = NA, nonpar = FALSE, thresh = 0.05, all_results = FALSE) 
# {
#     if (!nonpar) {
#         aov.nm <- "One-way ANOVA"
#     }
#     else {
#         aov.nm <- "Kruskal Wallis Test"
#     }
#     sig.num <- 0
#     aov.res <- sig.mat <- NULL
#     data <- as.matrix(mSetObj$dataSet$norm)
#     cls <- mSetObj$dataSet$cls
#     my.res <- GetFtestRes(mSetObj, nonpar)
#     aov.res <- my.res$aov.res
#     res <- my.res$f.res
#     fstat <- res[, 1]
#     p.value <- res[, 2]
#     names(fstat) <- names(p.value) <- colnames(data)
#     fdr.p <- p.adjust(p.value, "fdr")
#     all.mat <- data.frame(signif(fstat, 5), signif(p.value, 
#                                                    5), signif(-log10(p.value), 5), signif(fdr.p, 5))
#     rownames(all.mat) <- names(p.value)
#     colnames(all.mat) <- c("F.stat", "p.value", "-log10(p)", 
#                            "FDR")
#     ord.inx <- order(p.value)
#     ord.mat <- all.mat[ord.inx, ]
#     write.csv(ord.mat, "anova_all_results.csv")
#     inx.imp <- fdr.p <= thresh
#     sig.num <- sum(inx.imp, na.rm = TRUE)
#     cat("A total of", sig.num, "significant features were found.")
#     res <- 0
#     sig.f <- sig.p <- sig.fdr <- 1
#     if (sig.num > 0) {
#         res <- 1
#         sig.f <- fstat[inx.imp]
#         sig.p <- p.value[inx.imp]
#         sig.fdr <- fdr.p[inx.imp]
#         if (exists("aov.res")) {
#             qs::qsave(aov.res[inx.imp], file = "aov_res_imp.qs")
#         }
#         sig.mat <- all.mat[inx.imp, ]
#         ord.inx <- order(sig.mat$p.value)
#         sig.mat <- sig.mat[ord.inx, ]
#     }
#     else {
#         sig.mat <- ord.mat[1:10, ]
#     }
#     aov <- list(aov.nm = aov.nm, nonpar = nonpar, sig.num = sig.num, 
#                 raw.thresh = thresh, thresh = -log10(thresh), p.value = p.value, 
#                 p.log = -log10(p.value), fdr.p = fdr.p, inx.imp = inx.imp, 
#                 sig.f = sig.f, sig.p = sig.p, sig.fdr = sig.fdr, sig.mat = sig.mat,
#                 aov.res = aov.res)
#     mSetObj$analSet$aov <- aov
#     return(mSetObj)
# }
# 
# Do posthoc tests on significant features from ANOVA tests
Calculate.ANOVA.posthoc <- function(mSetObj=NA, post.hoc="fisher", thresh=0.05){
    
    sig.num <- mSetObj$analSet$aov$sig.num;
    inx.imp <- mSetObj$analSet$aov$inx.imp;
    sig.f <- mSetObj$analSet$aov$sig.f;
    sig.p <- mSetObj$analSet$aov$sig.p;
    sig.fdr <- mSetObj$analSet$aov$sig.fdr;
    nonpar <- mSetObj$analSet$aov$nonpar;
    cmp.res <- NULL;
    post.nm <- NULL;
    
    if(nonpar){
        sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(sig.fdr,5), 'NA');
        colnames(sig.mat) <- c("chi.squared", "p.value", "-log10(p)", "FDR", "Post-Hoc");
        fileName <- "kw_posthoc.csv";
    }else{     
        fileName <- "anova_posthoc.csv";    
        
        # do post-hoc only for signficant entries
        # note aov obj is not avaible using fast version
        # need to recompute using slower version for the sig ones
         aov.imp <- qs::qread("aov_res_imp.qs");

        # note this is only for post-hoc analysis. max 1000 in case too large
        if(sig.num > 1000){
            # update inx.imp   
            my.ranks <- rank(sig.p);
            inx.imp <- my.ranks <= 1000;
            aov.imp <- aov.imp[inx.imp];
        }
        
        if(post.hoc=="tukey"){
            tukey.res<-lapply(aov.imp, TukeyHSD, conf.level=1-thresh);
            my.cmp.res <- unlist(lapply(tukey.res, parseTukey, cut.off=thresh));
            post.nm = "Tukey's HSD";
        }else{
            fisher.res<-lapply(aov.imp, FisherLSD, thresh);
            my.cmp.res <- unlist(lapply(fisher.res, parseFisher, cut.off=thresh));
            post.nm = "Fisher's LSD";
        }
        
        cmp.res <- my.cmp.res;
        # post hoc only top 1000;
        
        if(sig.num > 1000){
            cmp.res <- rep(NA, sig.num); 
            cmp.res[inx.imp] <- my.cmp.res;
            post.nm <- paste(post.nm, "(top 1000)");
        }
        # create the result dataframe,
        # note, the last column is string, not double
        sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(sig.fdr,5), cmp.res);
        colnames(sig.mat) <- c("f.value", "p.value", "-log10(p)", "FDR", post.nm);
    }
    
    rownames(sig.mat) <- names(sig.p);
    # order the result simultaneously
    ord.inx <- order(sig.p, decreasing = FALSE);
    sig.mat <- sig.mat[ord.inx,,drop=F];
    
    # note only display top 1000 max for web (save all to the file)
    write.csv(sig.mat,file=fileName);
    if(sig.num > 1000){
        sig.mat <- sig.mat[1:1000,];
    }
    
    aov <- mSetObj$analSet$aov; 
    # add to the list, don't use append, as it does not overwrite
    aov$sig.nm <- fileName;
    aov$post.hoc <- post.hoc;
    aov$sig.mat <- sig.mat;
    
    mSetObj$analSet$aov <- aov;
    return(mSetObj);  
}
