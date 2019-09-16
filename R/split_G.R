#' Output the maximum potential scale reduction statistic of all parameters
#' estimated
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param Kfold an integer number of the number of cross-validation folds.
#'
#' @return a matrix made of perplexity from the results of cross-validation.
#'
#'
#' @importFrom magrittr %>%
#' @importFrom gtools mixedsort
#'
#' @examples
#' 
#' load(system.file("extdata/sample.rdata", package = "selectKSigs"))
#' G_split <- splitG(G, Kfold = 3)
#'
#' @export
#

splitG <- function(inputG, Kfold = 3) {
    # read in mutation feature data
    feature <- apply(inputG@featureVectorList, 2, paste, collapse="")
    count <- sum(inputG@countData[3,])

    # create a primary key for every single mutation
    f_s <- rep(paste0(rep(feature,
                          table(inputG@countData[1,])), "_",
                          inputG@countData[2,]),
               inputG@countData[3,])

    # randomly shuffle the order of mutations
    f_s <- f_s[sample(seq_len(length(f_s)))]
    # split mutations into Kfold number of folds
    folds <- cut(seq(1,length(f_s)),breaks=Kfold,labels=FALSE)

    # create a list to store results
    Glist <-  as.list(rep(NA, Kfold))

    # the splitting process is conducted for Kfold times
    for(i in seq_len(Kfold)){

      # then save the training and test data to the ith of Glist
      Glist[[i]] <- 
        list(select_kth_fold(inputG, k=i, f_s, folds, include = FALSE),
             select_kth_fold(inputG, k=i, f_s, folds, include = TRUE))
    }

    # return the entire Glist
    return(Glist)

}



