computeNoiseForSubset <- function (data, sample.size, n.iter) {
    # initialize matrices for storing intrinsic and extrinsic noise estimates
    int <- matrix (0, ncol=n.iter, nrow=7)
    ext <- matrix (0, ncol=n.iter, nrow=4)
    
    # in each iteration
    # sample a subset of cells
    # and compute the estimates
    for (i in 1:n.iter){
        subset.id <- sample (1:nrow(data), sample.size)
        int[,i] <- unlist (computeIntrinsicNoise (data[subset.id,1], data[subset.id,2])) / mean (data[subset.id,1]) / mean (data[subset.id,1])
        ext[,i] <- unlist (computeExtrinsicNoise (data[subset.id,1], data[subset.id,2])) / mean (data[subset.id,1]) / mean (data[subset.id,1])
    }
    
    rownames (int) <- c("ELSS", "unbiasedGeneral", "unbiasedEqualMean", "minMSEGeneral", "minMSEEqualMean", "asymptoticGeneral", "asymptoticEqualMean")
    rownames (ext) <- c("ELSS", "unbiased", "minMSE", "asymptotic")
    
    # compute mean and standard deviation
    # of estimates from multiple iterations
    int.mean <- apply (int, 1, mean)
    int.sd <- apply (int, 1, sd)
    ext.mean <- apply (ext, 1, mean)
    ext.sd <- apply (ext, 1, sd)
    
    # return estimates and summary statistics
    return (list(intrinsic=int, extrinsic=ext, intrinsic.mean=int.mean, intrinsic.sd=int.sd, extrinsic.mean=ext.mean, extrinsic.sd=ext.sd))
}

