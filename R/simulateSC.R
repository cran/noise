simulateSC <- function (n=1000, intrinsic=0.7, extrinsic=0.8, mean=1) {
    # simulate mean expression level for each cell
    # variance across cells is extrinsic noise
    mean.exp <- rnorm (n, mean=1, sd=sqrt(extrinsic))
    # simulate expression level of each reporter
    # conditioned on cell-specific mean and variance (intrinsic noise)
    reporter1.simu <- rnorm (n, mean=mean.exp, sd=sqrt(intrinsic))
    reporter2.simu <- rnorm (n, mean=mean.exp, sd=sqrt(intrinsic))
    
    return (data.frame (reporter1=reporter1.simu, reporter2=reporter2.simu))
}
