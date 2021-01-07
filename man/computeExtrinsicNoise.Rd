\name{computeExtrinsicNoise}
\alias{computeExtrinsicNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute estimates of the extrinsic noise.
}
\description{
This function computes several estimates of the extrinsic noise (unscaled by the mean).  The estimators, described in Fu and Pachter (2016), include the original estimators developed in Elowitz et al. (2002), the unbiased estimator, an min-MSE estimator, and an asymptotic estimator for large sample sizes.
}
\usage{
computeExtrinsicNoise(reporter1, reporter2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reporter1}{
A vector of continuous values.  
}
  \item{reporter2}{
A vector of continuous values.
}
}
\value{
Four (unscaled) estimates of extrinsic noise: the original estimators developed in Elowitz et al. (2002), the unbiased estimator, an min-MSE estimator, and an asymptotic estimator for large sample sizes.
}
\references{
Fu, A. Q. and Pachter, L. (2016). Estimating intrinsic and extrinsic noise from single-cell gene expression measurements.  arXiv:1601.03334.
Elowitz, M. B., A. J. Levine, E. D. Siggia, and P. S. Swain (2002) Stochastic gene expression in a single cell. Science, 297, 1183-1186.
}
\author{
Audrey Q. Fu
}
\seealso{
\code{\link{computeIntrinsicNoise}}, \code{\link{simulateSC}}.  See estimates for data \code{\link{elowitz_data}} and \code{\link{yang_nl10}}.
}

