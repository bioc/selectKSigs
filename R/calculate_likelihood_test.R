#' Output the maximum potential scale reduction statistic of all parameters
#' estimated
#'
#' @param train a MutationFeatureData S4 class output of training data.
#' @param test a MutationFeatureData S4 class output of test data.
#' @param param
#'
#'
#' @return the likelihood of the test data
#'
#'

Calculate_Likelihood_test <- function(train = G_train, test = G_test,
                                      param = param_train){

    F <- param@signatureFeatureDistribution
    Q <- param@sampleSignatureDistribution
    fdim <- slot(train, "possibleFeatures")
    K <- param@signatureNum
    isBG <- param@isBackGround
    BG <- 0
    sampleNum <- length(slot(train, "sampleList"))
    tol <- 1e-4

    p0 <- c(convertToTurbo_F(as.vector(F), fdim, K, isBG),
            convertToTurbo_Q(as.vector(Q), K, sampleNum))
    Y <- list(list(sampleNum, fdim, slot(test, "featureVectorList"),
                   slot(test, "countData")), K, isBG, BG)

    return(calcPMSLikelihood(p0, Y))
}

#' A function for calculating the log-likelihood from the data and parameters
#' 
#' @param p this variable includes the parameters for mutation signatures and 
#'          membership parameters
#' @param y this variable includes the information on the mutation features, 
#'          the number of mutation signatures specified and so on
#' @return a value
#' 
calcPMSLikelihood <- function(p, y) {
    
    sampleNum <- y[[1]][[1]]
    fdim <- y[[1]][[2]]
    patternList <- y[[1]][[3]]
    sparseCount <- y[[1]][[4]]
    K <- y[[2]]
    isBG <- y[[3]]
    BG0 <- y[[4]]
    
    patternNum <- ncol(patternList)
    samplePatternNum <- ncol(sparseCount)
    
    
    if (isBG) {
        varK <- K - 1
    } else {
        varK <- K
    }
    
    lenF <- varK * (sum(fdim) - length(fdim))
    lenQ <- (K - 1) * sampleNum
    F <- convertFromTurbo_F(p[seq_len(lenF)], fdim, K, isBG)
    Q <- convertFromTurbo_Q(p[(lenF + 1):(lenF + lenQ)], K, sampleNum)
    
    dim(Q) <- c(sampleNum, K)
    Q <- t(Q)
    ####################
    
    return(getLogLikelihoodC(as.vector(patternList), as.vector(sparseCount), 
                             as.vector(F), as.vector(Q), fdim, K, 
                             sampleNum, patternNum, samplePatternNum, 
                             isBG, BG0))
    
}
