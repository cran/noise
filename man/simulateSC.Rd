\name{simulateSC}
\alias{simulateSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate single-cell expression levels of two reporters
}
\description{This function simulates expression levels of two reporters across 
single cells, mimicking the two-reporter assay.  The hierarchical model described in Fu 
and Pachter (2016) is used for simulation.  We further make the simplifying assumption 
that intrinsic noise is the same across cells.
}
\usage{
simulateSC(n = 1000, intrinsic = 0.7, extrinsic = 0.8, mean = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of single cells (sample size).
}
  \item{intrinsic}{Scalar.  The (unscaled) intrinsic noise (or within-cell variability), denoted 
by \eqn{\sigma^2} (equation 8) in Fu and Pachter (2016).
}
  \item{extrinsic}{Scalar.  The (unscaled) extrinsic noise (or between-cell variability), denoted 
by \eqn{\sigma^2_\mu} (equation 9) in Fu and Pachter (2016).
}
  \item{mean}{Scalar.  The overall mean of expression level, denoted by \eqn{\mu} (equation 6) in Fu and Pachter (2016).
}
}
\value{
A data frame of two columns and \eqn{n} rows.  Each column contains the expression levels of a reporter.  Each row is a single cell.
}
\references{
Fu, A. Q. and Pachter, L. (2016). Estimating intrinsic and extrinsic noise from single-cell gene expression measurements.  arXiv:1601.03334.
}
\author{
Audrey Q. Fu
}
\seealso{
\code{\link{computeIntrinsicNoise}}, \code{\link{computeExtrinsicNoise}}.
}
\examples{
# simulation #1

# simulate 500 data sets
n.simu <- 500

# true intrinsic and extrinsic noise
int.true <- 0.7
ext.true <- 0.8

# create matrices to hold estimated intrinsic and extrinsic noise
# using different estimators
int.simu.mtx <- matrix (0, nrow=n.simu, ncol=8)
ext.simu.mtx <- matrix (0, nrow=n.simu, ncol=4)

for (i in 1:n.simu) {
    n <- 1000
    simu <- simulateSC (n=n, intrinsic=int.true, extrinsic=ext.true, mean=1)
    
    int.simu.mtx[i,] <- c(unlist (computeIntrinsicNoise (simu[,1], simu[,2])), 
        cor (simu[,1], simu[,2]))
    ext.simu.mtx[i,] <- unlist (computeExtrinsicNoise (simu[,1], simu[,2]))
    
}

# add column names to simulation estimates
colnames (int.simu.mtx) <- c("ELSS", "unbiasedGeneral", "unbiasedEqualMean", 
    "minMSEGeneral", "minMSEEqualMean", "asymptoticGeneral", 
    "asymptoticEqualMean", "cor")
colnames (ext.simu.mtx) <- c("ELSS", "unbiased", "minMSE", "asymptotic")


# simulation #2

# simulate 500 data sets
n.simu <- 500

# true intrinsic and extrinsic noise
int.true <- 0.7
ext.true <- 0.8

# use true correlation for the min-MSE estimates of extrinsic noise
true.cor <- ext.true / (ext.true + int.true)

# create matrices to hold estimated intrinsic and extrinsic noise
# using different estimators
int.simu.mtx <- matrix (0, nrow=n.simu, ncol=8)
ext.simu.mtx <- matrix (0, nrow=n.simu, ncol=4)
ext.simu.mtx.2 <- matrix (0, nrow=n.simu, ncol=4)

for (i in 1:n.simu) {
    n <- 50
    simu <- simulateSC (n=n, intrinsic=int.true, extrinsic=ext.true, mean=1)
    
    int.simu.mtx[i,] <- c(unlist (computeIntrinsicNoise (simu[,1], simu[,2])), 
      cor (simu[,1], simu[,2]))
    ext.simu.mtx[i,] <- unlist (computeExtrinsicNoise (simu[,1], simu[,2]))
    ext.simu.mtx.2[i,] <- c(unlist (computeExtrinsicNoiseKnownCor (simu[,1], 
      simu[,2], true.cor)))
}

# add column names to simulation estimates
colnames (int.simu.mtx) <- c("ELSS", "unbiasedGeneral", "unbiasedEqualMean", 
    "minMSEGeneral", "minMSEEqualMean", "asymptoticGeneral", 
    "asymptoticEqualMean", "cor")
colnames (ext.simu.mtx) <- c("ELSS", "unbiased", "minMSE", "asymptotic")
colnames (ext.simu.mtx.2) <- c("ELSS", "unbiased", "minMSE", "asymptotic")

# compute the MSE of estimates
computeMSE <- function (a, t) {return (mean((a-t)^2))}
apply (int.simu.mtx[,1:7], 2, computeMSE, t=int.true)
apply (ext.simu.mtx, 2, computeMSE, t=ext.true)
apply (ext.simu.mtx.2, 2, computeMSE, t=ext.true)

}

