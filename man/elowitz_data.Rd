\name{elowitz_data}
\alias{elowitz_data}
\docType{data}
\title{
Expression of reporter genes in the D22 and M22 E. coli cells from Elowitz et al (2002).
}
\description{
Expression of reporter genes CFP and YFP in over 200 E. coli cells of two strains: D22 and M22.  These values are displayed in a scatterplot in Elowitz et al (2002) Fig 3a.
}
\usage{data("elowitz_data")}
\format{
  The format is:
List of 2
 $ D22:'data.frame':	284 obs. of  2 variables:
  ..$ CFP: num [1:284] 3080 3082 2893 3053 2891 ...
  ..$ YFP: num [1:284] 2309 2394 2145 2340 2245 ...
 $ M22:'data.frame':	250 obs. of  2 variables:
  ..$ CFP: num [1:250] 2438 2316 2521 2646 2830 ...
  ..$ YFP: num [1:250] 1409 1391 1511 1460 1638 ...
}
\references{
Elowitz, M. B., A. J. Levine, E. D. Siggia, and P. S. Swain (2002) Stochastic gene expression in a single cell. Science, 297, 1183-1186.
}
\examples{
data(elowitz_data)

# Normalize data such that they are 
# comparable to Fig 3a in Elowitz et al. (2002).
# Normalized data have mean 1.
D22.cfp.norm <- (elowitz_data$D22[,1]-mean (elowitz_data$D22[,1]))/sd(elowitz_data$D22[,1])/8+1
D22.yfp.norm <- (elowitz_data$D22[,2]-mean (elowitz_data$D22[,2]))/sd(elowitz_data$D22[,2])/8+1

M22.cfp.norm <- (elowitz_data$M22[,1]-mean (elowitz_data$M22[,1]))/sd(elowitz_data$M22[,1])/12+1
M22.yfp.norm <- (elowitz_data$M22[,2]-mean (elowitz_data$M22[,2]))/sd(elowitz_data$M22[,2])/12+1

# Compute noise estimates.
# Since the mean is 1, estimates with and without 
# the scaling are the same.
unlist (computeIntrinsicNoise (D22.cfp.norm, D22.yfp.norm))
unlist (computeExtrinsicNoise (D22.cfp.norm, D22.yfp.norm))

unlist (computeIntrinsicNoise (M22.cfp.norm, M22.yfp.norm))
unlist (computeExtrinsicNoise (M22.cfp.norm, M22.yfp.norm))

}
\keyword{datasets}
