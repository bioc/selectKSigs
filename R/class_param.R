#' Get the values of loglikelihood
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return likelihood values estimated by \code{pmgetSignature} in \code{HiLDA}
#' 
#' 
setGeneric("getLL", function(object) standardGeneric("getLL"))
setMethod("getLL", "EstimatedParameters", function(object) object@loglikelihood)


#' Get the number of signatures
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return the number of signatures in \code{pmgetSignature} in \code{HiLDA}
#' 
#' 
setGeneric("getK", function(object) standardGeneric("getK"))
setMethod("getK", "EstimatedParameters", function(object) object@signatureNum)


#' Get the statsus of using the background signature
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return the status of using the background signature
#' 
#' 
setGeneric("getBG", function(object) standardGeneric("getBG"))
setMethod("getBG", "EstimatedParameters", function(object) object@isBackGround)


#' Get the sample list
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return the sample list of named elements. 
#' 
#' 
setGeneric("getSamplelist", function(object) standardGeneric("getSamplelist"))
setMethod("getSamplelist", "EstimatedParameters", 
          function(object) object@sampleList)


#' Get a vector of possible features
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return a vector of possible features
#' 
#' 
setGeneric("getFeatures", function(object) standardGeneric("getFeatures"))
setMethod("getFeatures", "EstimatedParameters", 
          function(object) object@possibleFeatures)


#' Get an array of signature feature distributions
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return an array of signature feature distributions
#' 
#' 
setGeneric("getSignatures", function(object) standardGeneric("getSignatures"))
setMethod("getSignatures", "EstimatedParameters", 
          function(object) object@signatureFeatureDistribution)


#' Get a matrix of mutational exposures of signatures
#' 
#' @param object the EstimatedParameters class (the result of pmgetSignature)
#' 
#' @return a matrix of mutational exposures of signatures
#' 
#' 
setGeneric("getExposures", function(object) standardGeneric("getExposures"))
setMethod("getExposures", "EstimatedParameters", 
          function(object) object@sampleSignatureDistribution)


