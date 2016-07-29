computeExtrinsicNoiseKnownCor <- function (reporter1, reporter2, true.cor) {
  # input: reporter1, reporter2: vectors of normalized reporter gene expression
  
  # compute sample size
  n <- length (reporter1)
  # compute \sum C_i Y_i - n\bar{C}\bar{Y}
  sum.cov <- sum (reporter1*reporter2) - n * mean (reporter1) * mean (reporter2)
  # compute sample correlation
  #sample.cor <- cor (reporter1, reporter2)
  
  # unbiased estimator
  unbiased <- sum.cov / (n-1)
  
  # min MSE estimator
  a <- 1/(true.cor^2) + (n-1)*(1+1/n)
  min.mse <- sum.cov / a
  
  # asymptotic estimator
  asym <- sum.cov / n
  
  # return a list of estimators
  return (list (ELSS=asym, unbiased=unbiased, minMSE=min.mse, asymptotic=asym))
}
