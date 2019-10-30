#' Get the statsus of specifying the transcription bias
#' 
#' @param object the MutationFeatureData class
#' 
#' @return the status of specifying the transcription bias
#' 
#' @examples print(getTranscription(G))
#' 
setGeneric("getTranscription", function(object) 
  standardGeneric("getTranscription"))
setMethod("getTranscription", "MutationFeatureData", 
          function(object) object@transcriptionDirection)


#' Get the sample list
#' 
#' @param object the MutationFeatureData class
#' 
#' @return the sample list of named elements. 
#' 
#' @examples print(getSamplelist(G))
#' 
setGeneric("getSamplelist", function(object) standardGeneric("getSamplelist"))
setGeneric("getSamplelist<-", 
           function(x, value) standardGeneric("getSamplelist<-"))

setMethod("getSamplelist", "MutationFeatureData", 
          function(object) object@sampleList)
setMethod("getSamplelist<-", "MutationFeatureData", function(x, value) {
  x@sampleList <- value
  x
})


#' Get the count data in a matrix
#' 
#' @param object the MutationFeatureData class
#' 
#' @return the count data in a matrix
#' 
#' @examples print(getCounts(G))
#' 
setGeneric("getCounts", function(object) standardGeneric("getCounts"))
setGeneric("getCounts<-", 
           function(x, value) standardGeneric("getCounts<-"))

setMethod("getCounts", "MutationFeatureData", function(object) object@countData)
setMethod("getCounts<-", "MutationFeatureData", function(x, value) {
  x@countData <- value
  x
})

#' Get a matrix of feature vector list
#' 
#' @param object the MutationFeatureData class
#' 
#' @return a matrix of feature vector list
#' 
#' @examples print(getFeatureVec(G))
#' 
setGeneric("getFeatureVec", function(object) standardGeneric("getFeatureVec"))
setGeneric("getFeatureVec<-", 
           function(x, value) standardGeneric("getFeatureVec<-"))


setMethod("getFeatureVec", "MutationFeatureData", 
          function(object) object@featureVectorList)
setMethod("getFeatureVec<-", "MutationFeatureData", function(x, value) {
  x@featureVectorList <- value
  x
})
