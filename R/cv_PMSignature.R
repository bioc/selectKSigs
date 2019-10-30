#' Output the maximum potential scale reduction statistic of all parameters
#' estimated
#'
#' @param inputG a MutationFeatureData S4 class.
#' @param Kfold an integer number of the number of cross-validation folds.
#' @param nRep an integer number of replications. 
#' @param Klimit an integer of the maximum value of number of signatures.
#'
#' @importFrom HiLDA pmgetSignature
#' 
#' @examples
#'
#' load(system.file("extdata/sample.rdata", package = "selectKSigs"))
#' results <- cv_PMSignature(G, Kfold = 3)
#' 
#' @return a matrix of measures
#'
#' @export
#'

cv_PMSignature <- function(inputG, Kfold = 3, nRep = 3, Klimit = 8){
    # create a list to store the likelihood
    LL_list <- replicate(nRep, matrix(NA, Klimit-1, Kfold), simplify=FALSE)
    # loop over the number of replications
    for(r in seq_len(nRep)){
      G_split <- splitG(inputG,  Kfold = Kfold)

      # repeat over a range of K values
      for(k in seq_len(Klimit-1)+1){
        # repeat over the folds
        for(i in seq_len(Kfold)){
          G_train <- G_split[[i]][[1]]
          G_test <- G_split[[i]][[2]]
          param_train <- HiLDA::pmgetSignature(G_train, K = k)

          LL_list[[r]][k-1, i] <- Calculate_Likelihood_test(G_train, G_test,
                                                            param_train)
        }
      }
    }
    # create a matrix to store the measures
    measuremat <- data.frame(matrix(NA, Klimit-1, 4))
    colnames(measuremat) <- c("AIC", "AICc", "BIC", "Perplexity")

    for(k in seq_len(Klimit-1)+1){
      param<-HiLDA::pmgetSignature(inputG, K = k)

      # number of parameters
      nP <- length(getSamplelist(param))*(k-1) +
        k*(sum(getFeatures(param))-length(getFeatures(param)))

      # AIC
      measuremat[k-1, 1]<- getLL(param)*(-2) + 2*nP

      #AICc
      measuremat[k-1, 2] <- getLL(param)*(-2) +
        2*nP + 2*(nP)*(nP+1)/(sum(getCounts(inputG)[3,]) - nP-1)

      # BIC
      measuremat[k-1, 3] <- getLL(param)*(-2) +
        log(sum(getCounts(inputG)[3,]))*nP

    }
      # Perplexity
      measuremat[, 4] <- exp(rowMeans(sapply(LL_list, function(x) rowMeans(x)))*
                               (-1)/sum(getCounts(inputG)[3,]))
    
    return(measuremat)
}
