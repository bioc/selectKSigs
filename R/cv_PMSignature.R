#' Output the maximum potential scale reduction statistic of all parameters
#' estimated
#'
#' @param train a MutationFeatureData S4 class output of training data.
#' @param test a MutationFeatureData S4 class output of test data.
#' @param param
#'
#' @importFrom HiLDA pmgetSignature
#' 
#' @examples
#'
#' res <- cv_PMSignature(G, Kfold = 3)
#' 
#' @return the likelihood of the test data
#'
#'

cv_PMSignature <- function(inputG, Kfold = 3, nRep = 3, Klimit = 8){
    # create a list to store the likelihood
    LL_list <- replicate(nRep, matrix(NA, Klimit-1, Kfold), simplify=F)

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
    measuremat <- matrix(NA, Klimit-1, 4)

    for(k in seq_len(Klimit-1)+1){
      param<-HiLDA::pmgetSignature(G,K = k)

      # number of parameters
      nP <- length(param@sampleList)*(k-1) +
        k*(sum(param@possibleFeatures)-length(param@possibleFeatures))

      # AIC
      measuremat[k-1, 1]<- param@loglikelihood*(-2) + 2*nP

      #AICc
      measuremat[k-1, 2] <-param@loglikelihood*(-2) +
        2*nP + 2*(nP)*(nP+1)/(sum(G@countData[3,]) - nP-1)

      # BIC
      measuremat[k-1, 3] <- param@loglikelihood*(-2) +
        log(sum(G@countData[3,]))*nP

      # Perplexity
      measuremat[k-1, 4] <- exp(param@loglikelihood*(-1)/sum(G@countData[3,]))
    }



}
