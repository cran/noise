\name{computeNoiseForSubset}
\alias{computeNoiseForSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select a random subset of data for noise estimation.
}
\description{
This function randomly selects a subset of cells (rows) from the data set, computes multiple estimates of intrinsic and extrinsic noise, as well as their mean and standard deviation.
}
\usage{
computeNoiseForSubset(data, sample.size, n.iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A numeric matrix of two columns.  Each row is a cell, and each column expression of a reporter gene.
}
  \item{sample.size}{
An integer that specifies the number of cells in the subset.
}
  \item{n.iter}{
An integer that specifies the number of iterations (for calcuation of mean and standard deviation).
}
}
\value{
A list that consists of the following components:
  \item{intrinsic}{A numeric matrix of esimated intrinsic noise. 7 rows and n.iter columns.}
  \item{extrinsic}{A numeric matrix of esimated extrinsic noise. 4 rows and n.iter columns.}
  \item{intrinsic.mean}{A numeric vector of length 7 that contains the mean estimates of intrinsic noise.}
  \item{intrinsic.sd}{A numeric vector of length 7 that contains the standard deviation of the estimates of intrinsic noise.}
  \item{extrinsic.mean}{A numeric vector of length 7 that contains the mean estimates of extrinsic noise.}
  \item{extrinsic.sd}{A numeric vector of length 7 that contains the standard deviation of the estimates of extrinsic noise.}
}
\references{
Fu, A. Q. and Pachter, L. (2016). Estimating intrinsic and extrinsic noise from single-cell gene expression measurements.  arXiv:1601.03334.
}
\author{
Audrey Q. Fu
}
\seealso{
\code{\link{computeIntrinsicNoise}}, \code{\link{computeExtrinsicNoise}}, \code{\link{elowitz_data}}, \code{\link{yang_nl10}}.
}
\examples{
data(yang_nl10)

# quantile normalization on log2 transformed data
# install bioconductor package for quantile normalization
# source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore')
library(preprocessCore)

# ignore a few values that are negative
yang_nl10.pos <- yang_nl10[-which (yang_nl10[,1]<0),]
yang_nl10.pos.log2.quant <- normalize.quantiles (as.matrix (log2 (yang_nl10.pos[,c(1,3)])))

# subset the data and compute noise estimates
yang.50 <- computeNoiseForSubset (yang_nl10.pos.log2.quant, sample.size=50, n.iter=1000)
summary (yang.50)
}