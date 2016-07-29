\name{computeIntrinsicNoise}
\alias{computeIntrinsicNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute estimates of the intrinsic noise.
}
\description{
This function computes several estimates of the intrinsic noise (unscaled by the mean).  The estimators, described in Fu and Pachter (2016), include the original estimators developed in Elowitz et al. (2002), unbiased estimators with and without assuming equal mean of the two reporters, min-MSE estimators with and without assuming equal mean, and asymptotic estimators for large sample sizes with and without assuming equal mean.
}
\usage{
computeIntrinsicNoise(reporter1, reporter2)
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
Six (unscaled) estimates of intrinsic noise: the original estimators developed in Elowitz et al. (2002), unbiased estimators with and without assuming equal mean of the two reporters, min-MSE estimators with and without assuming equal mean, and asymptotic estimators for large sample sizes with and without assuming equal mean.
}
\references{
Fu, A. Q. and Pachter, L. (2016). Estimating intrinsic and extrinsic noise from single-cell gene expression measurements.  arXiv:1601.03334.
Elowitz, M. B., A. J. Levine, E. D. Siggia, and P. S. Swain (2002) Stochastic gene expression in a single cell. Science, 297, 1183-1186.
}
\author{
Audrey Q. Fu
}
\seealso{
\code{\link{computeExtrinsicNoise}}, \code{\link{simulateSC}}.  See estimates for data \code{\link{elowitz_data}} and \code{\link{yang_nl10}}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
