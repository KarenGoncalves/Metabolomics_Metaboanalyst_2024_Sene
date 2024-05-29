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