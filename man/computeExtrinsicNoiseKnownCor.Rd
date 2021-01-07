\name{computeExtrinsicNoiseKnownCor}
\alias{computeExtrinsicNoiseKnownCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute estimates of the extrinsic noise using a known correlation in the min-MSE estimator.
}
\description{
This function is similar to \code{\link{computeExtrinsicNoise}}, and computes several estimates of the extrinsic noise (unscaled by the mean).  The estimators, described in Fu and Pachter (2016), include the original estimators developed in Elowitz et al. (2002), the unbiased estimator, an min-MSE estimator, and an asymptotic estimator for large sample sizes.  The only difference between this function calculates the min-MSE estimate using a given correlation.
}
\usage{
computeExtrinsicNoiseKnownCor(reporter1, reporter2, true.cor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reporter1}{
A vector of continuous values.  
}
  \item{reporter2}{
A vector of continuous values.
}
  \item{true.cor}{
A scalar.
}
}
\value{
Four (unscaled) estimates of extrinsic noise: the original estimators developed in Elowitz et al. (2002), the unbiased estimator, an min-MSE estimator (using the given correlation), and an asymptotic estimator for large sample sizes.
}
\references{
Fu, A. Q. and Pachter, L. (2016). Estimating intrinsic and extrinsic noise from single-cell gene expression measurements.  arXiv:1601.03334.
Elowitz, M. B., A. J. Levine, E. D. Siggia, and P. S. Swain (2002) Stochastic gene expression in a single cell. Science, 297, 1183-1186.
}
\author{
Audrey Q. Fu
}
\seealso{
\code{\link{computeExtrinsicNoise}}, \code{\link{simulateSC}}.
}
