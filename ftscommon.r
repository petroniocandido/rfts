FLR <- function(a,b){
    nc <- list(
        lhs = a,
        rhs = b
    )
    nc$dump <- function() return(paste(nc$lhs,nc$rhs,sep="->"))
    return (nc)
}


fuzzyInstance <- function(inst, fuzzySets){
    mv <- list()
    for(set in fuzzySets) mv[[set[["name"]]]] <- set$membership(inst) 
    return (mv)
}


fuzzySeries <- function(data,fuzzySets){
    fuzzydata <- c()
    for(i in 1:length(data)){
        mv <- fuzzyInstance(data[i],fuzzySets)
        mvals <- unlist(mv, use.names=FALSE)
        best_set <- which(mvals == max(mvals))[1]
        fuzzydata[i] <- names(mv)[best_set]
    }
    return (fuzzydata)
}


genFLR <- function(fuzzyData){
    flrs <- list();
    for(i in 2:length(fuzzyData)){
        tmp <- FLR(fuzzyData[i-1],fuzzyData[i]);
        if(is.null(flrs[[tmp$dump()]])) flrs[[tmp$dump()]] <- tmp;
    }
    return (flrs);
}

genRecurrentFLR <- function(fuzzyData){
    flrs <- list();
    for(i in 2:length(fuzzyData)){
        tmp <- FLR(fuzzyData[i-1],fuzzyData[i]);
        flrs[[i-1]] <- tmp;
    }
    return (flrs);
}

forecastAhead <- function(fts, x,n){
        forecasts <- rep(0,n)
        forecasts[1] <- fts$forecast(x)
        for(i in 2:n) {
			forecasts[i] <- fts$forecast(forecasts[i-1])
		}
		return (forecasts)
 }
