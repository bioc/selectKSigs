#include <Rcpp.h>
using namespace Rcpp;


//' Calculate the value of the log-likelihood for given parameters
//'
//' @param vPatternList The list of possible mutation features (converted
//' to a vector)
//' @param vSparseCount The table showing (mutation feature, sample, the number
//' of mutation) (converted to a vector)
//' @param vF F (converted to a vector)
//' @param vQ Q (converted to a vector)
//' @param fdim a vector specifying the number of possible values for each
//' mutation signature
//' @param signatureNum the number of mutation signatures
//' @param sampleNum the number of cancer genomes
//' @param patternNum the number of possible combinations of all the mutation
//' features
//' @param samplePatternNum the number of possible combination of samples and
//' mutation patternns
//' @param isBackground the logical value showing whether a background mutaiton
//' features is included or not
//' @param vF0 a background mutaiton features
//' @return a value
// [[Rcpp::export]]
double getLogLikelihoodC(NumericVector vPatternList, NumericVector vSparseCount,
                         NumericVector vF, NumericVector vQ, NumericVector fdim,
                         int signatureNum, int sampleNum, int patternNum,
                         int samplePatternNum, bool isBackground,
                         NumericVector vF0) {

  NumericVector vTheta(signatureNum * samplePatternNum);
  NumericVector vF_full(signatureNum * patternNum);


  int variableSigNum;
  if (isBackground) {
    variableSigNum = signatureNum - 1;
    for (int m = 0; m < patternNum; m++) {
      vF_full[signatureNum - 1 + m * signatureNum] = vF0[m];
    }
  } else {
    variableSigNum = signatureNum;
  }


  for (int m = 0; m < patternNum; m++) {
    for (int k = 0; k < variableSigNum; k++) {
      vF_full[k + m * signatureNum] = 1;
    }
  }


  int featureInd;
  for (int m = 0; m < patternNum; m++) {

    for (int l = 0; l < fdim.size(); l++) {
      featureInd = vPatternList[l + m * fdim.size()] - 1;
      for (int k = 0; k < variableSigNum; k++) {
        vF_full[k + m * signatureNum] =
          vF_full[k + m * signatureNum] * vF[k + l * variableSigNum +
          featureInd * variableSigNum * fdim.size()];
      }
    }

  }


  double tempSum = 0;
  int sampleInd;
  int patternInd;

  double logLikelihood = 0;
  for (int nm = 0; nm < samplePatternNum; nm++) {

    patternInd = vSparseCount[0 + nm * 3] - 1;
    sampleInd = vSparseCount[1 + nm * 3] - 1;

    tempSum = 0;
    for(int k = 0; k < signatureNum; k++) {
      tempSum = tempSum + vQ[k +
        sampleInd * signatureNum] * vF_full[k + patternInd * signatureNum];

    }

    if (tempSum > 1e-10 and vSparseCount[2 + nm * 3] > 0) {
      logLikelihood = logLikelihood + log(tempSum) * vSparseCount[2 + nm * 3];
    }

  }

  return(logLikelihood);

}




