#' Get the statsus of specifying the transcription bias
#' 
#' @param object the MutationFeatureData class
#' 
#' @return the status of specifying the transcription bias
#' 
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
#' 
setGeneric("getSamplelistG", function(object) standardGeneric("getSamplelistG"))
setGeneric("getSamplelistG<-", 
           function(x, value) standardGeneric("getSamplelistG<-"))

setMethod("getSamplelistG", "MutationFeatureData", 
          function(object) object@sampleList)
setMethod("getSamplelistG<-", "MutationFeatureData", function(x, value) {
  x@sampleList <- value
  x
})


#' Get the count data in a matrix
#' 
#' @param object the MutationFeatureData class
#' 
#' @return the count data in a matrix
#' 
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
