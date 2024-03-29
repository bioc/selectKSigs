% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_PMSignature.R
\name{cv_PMSignature}
\alias{cv_PMSignature}
\title{Output the maximum potential scale reduction statistic of all parameters
estimated}
\usage{
cv_PMSignature(inputG, Kfold = 3, nRep = 3, Klimit = 8)
}
\arguments{
\item{inputG}{a MutationFeatureData S4 class.}

\item{Kfold}{an integer number of the number of cross-validation folds.}

\item{nRep}{an integer number of replications.}

\item{Klimit}{an integer of the maximum value of number of signatures.}
}
\value{
a matrix of measures
}
\description{
Output the maximum potential scale reduction statistic of all parameters
estimated
}
\examples{

load(system.file("extdata/sample.rdata", package = "selectKSigs"))
results <- cv_PMSignature(G, Kfold = 3)

}
