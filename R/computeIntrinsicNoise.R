computeIntrinsicNoise <- function (reporter1, reporter2) {
    # input: reporter1, reporter2: vectors of normalized reporter gene expression
    
    # compute sample size
    n <- length (reporter1)
    # compute sum of squared differences
    sum.diff <- sum ((reporter1-reporter2)^2)
    # compute n*(\bar{C} - \bar{Y})^2
    sum.mean.diff <- n*((mean (reporter1) - mean (reporter2))^2)
    # compute sample correlation
    sample.cor <- cor (reporter1, reporter2)
    # compute estimated ratio of extrinsic vs intrinsic noise
    est.ratio <- sample.cor / (1-sample.cor)
    
    # unbiased estimator under general conditions
    unbiased.general <- (sum.diff - sum.mean.diff) / 2 / (n-1)
    # unbiased estimator assuming that the two repoters have equal means
    unbiased.eq.mean <- sum.diff / 2 / n
    
    # min MSE estimator under general conditions
    a <- (2*n^3 - 7*n + 6) / (2*(n^2 - n)) + (2 - n) / (n^2 - n) * est.ratio + 1/(2*(n^2 - n)) * est.ratio^2
    min.mse.general <- (sum.diff - sum.mean.diff) / 2 / a
    # min MSE estimator assuming that the two reporters have equal means
    min.mse.eq.mean <- sum.diff / 2 / (n+2)
    
    # asymptotic estimator
    asym.general <- (sum.diff - sum.mean.diff) / 2 / n
    
    # return a list of estimators
    return (list (ELSS=unbiased.eq.mean, unbiasedGeneral=unbiased.general, unbiasedEqualMean=unbiased.eq.mean, minMSEGeneral=min.mse.general, minMSEEqualMean=min.mse.eq.mean, asymptoticGeneral=asym.general, asymptoticEqualMean=unbiased.eq.mean))
}
