#' Output the training data or test data
#'
#' @param inputG a MutationFeatureData S4 class output by the pmsignature.
#' @param k an integer number of the number of cross-validation folds.
#' @param f_s a primary key of combining the feature pattern and sample ID.
#' @param folds the assignment to each fold.
#' @param include a boolean indictor of whether to include kth fold or not.
#'
#' @return a MutationFeatureData S4 class of either include or exclude kth fold.
#'
#' @importFrom magrittr %>%
#' @importFrom gtools mixedsort
#'
#


select_kth_fold <- function(inputG, k, f_s, folds, include){
    G_temp <- inputG
    # select or unselect the ith fold of inputG
    if (include == TRUE){
      # create test data
      temp <- f_s[which(folds==k,arr.ind=TRUE)] %>% sort()
    } else {
      # create traning data
      temp <- f_s[which(folds!=k,arr.ind=TRUE)] %>% sort()
    }
    # split the primary key to extract ID number and pattern
    temp_unlist <- strsplit(temp, "_")
    # t_list stores the feature pattern
    t_list <- as.numeric(lapply(temp_unlist, "[[", 1) %>% unique())
    # tabulate each mutation pattern, i.e., "124242"
    rawCount <- matrix(as.numeric(unlist(temp_unlist)), ncol=2, byrow=TRUE) %>%
      data.frame()
    rawCount[,c(2,1)] <- rawCount
    # return
    rawCount[,2] <- unlist(lapply(rawCount[,2], function(x) which(x == t_list)))

    # tabulate each mutation type
    tableCount <- table(rawCount)
    # only retain those are greater than zero
    w <- which(tableCount > 0, arr.ind=TRUE)
    # create a table with mutation feature, sample, counts
    procCount <- cbind(w[,2], w[,1], tableCount[w])
    rownames(procCount) <- NULL

    # determine whether its has transcription direction
    ncol = ifelse(getTranscription(inputG), 6, 5)

    # split the feature into a vector, i.e.,"124242" to "1, 2, 4, 2, 4, 2"
    getFeatureVec(G_temp) <-
      matrix(unlist(lapply(lapply(temp_unlist, "[[", 1) %>% unique(),
                           function(x) strsplit(x, ""))) %>%
               as.numeric(), ncol=ncol , byrow = TRUE) %>% t()
    getSamplelistG(G_temp) <-
      paste0("sample","_",unique(unlist(lapply(temp_unlist, "[[", 2))) %>%
               gtools::mixedsort())
    getCounts(G_temp)  <- t(procCount)

    return(G_temp)
}
