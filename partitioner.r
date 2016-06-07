UniversePartitioner <- function(data, npart, mf, prefix){
    fzy_sets <- list()
    bounds <- (max(data) - min(data))/npart
    dmax <- max(data) + bounds
    dmin <- min(data) - bounds
    dlen <- dmax - dmin
    partlen <- dlen / npart
    partition <- dmin
    for(c in 1:npart){
        setname <- paste(prefix,as.character(c),sep="")
        tmp <-FuzzySet(setname,mf,c(partition-partlen, partition, partition+partlen), partition )
        fzy_sets[[setname]] <-tmp 
        partition <- partition + partlen
    }
    return (fzy_sets)
}
