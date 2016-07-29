\name{yang_nl10}
\alias{yang_nl10}
\docType{data}
\title{
Expression of reporter genes in the NL010 E. coli cells from Yang et al (2014).
}
\description{
Expression of reporter genes CFP and mCherry in over 40,000 E. coli cells.  A subset of these values are displayed in a scatterplot in Yang et al (2014) Fig 3a rightmost panel.
}
\usage{data("yang_nl10")}
\format{
  A data frame with 40683 observations on the following 3 variables.
  \describe{
    \item{\code{CFP}}{a numeric vector}
    \item{\code{Venus}}{a numeric vector}
    \item{\code{mCherry}}{a numeric vector}
  }
}
\references{
Yang, S., S. Kim, Y. R. Lim, C. Kim, H. J. An, J.-H. Kim, J. Sung, and N. K. Lee (2014) Contribution of RNA polymerase concentration variation to protein expression noise. Nature Communications, 5, 4761.
}
\examples{
data(yang_nl10)

# compute the noise estimates

# no normalization
# unscaled by mean
unlist (computeIntrinsicNoise (yang_nl10[,1], yang_nl10[,3]))
unlist (computeExtrinsicNoise (yang_nl10[,1], yang_nl10[,3]))

# scaled by mean
unlist (computeIntrinsicNoise (yang_nl10[,1], yang_nl10[,3])) / 
    mean (yang_nl10[,1]) / mean(yang_nl10[,3])
unlist (computeExtrinsicNoise (yang_nl10[,1], yang_nl10[,3])) / 
    mean (yang_nl10[,1]) / mean(yang_nl10[,3])

# quantile normalization on log2 transformed data
# install bioconductor package for quantile normalization
# source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore')
library(preprocessCore)

# ignore a few values that are negative
yang_nl10.pos <- yang_nl10[-which (yang_nl10[,1]<0),]
yang_nl10.pos.log2.quant <- normalize.quantiles (as.matrix (log2 (yang_nl10.pos[,c(1,3)])))

# unscaled by mean
unlist (computeIntrinsicNoise (yang_nl10.pos.log2.quant[,1], yang_nl10.pos.log2.quant[,2]))
unlist (computeExtrinsicNoise (yang_nl10.pos.log2.quant[,1], yang_nl10.pos.log2.quant[,2]))

# scaled by mean
unlist (computeIntrinsicNoise (yang_nl10.pos.log2.quant[,1], yang_nl10.pos.log2.quant[,2])) / 
    mean (yang_nl10.pos.log2.quant[,1]) / mean(yang_nl10.pos.log2.quant[,2])
unlist (computeExtrinsicNoise (yang_nl10.pos.log2.quant[,1], yang_nl10.pos.log2.quant[,2])) / 
    mean (yang_nl10.pos.log2.quant[,1]) / mean(yang_nl10.pos.log2.quant[,2])

}
\keyword{datasets}
