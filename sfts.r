SFLRG <- function(seaz, prhs){
    nc <- list(
        season = seaz,
        rhs = prhs
    );
    
    nc$put <- function(x){
        if(!(x %in% nc$rhs)){
            return (SFLRG(nc$season, c(nc$rhs,x)));            
        } else { 
            return (nc);
        }
    } 
            
    nc$dump <- function() {
        prhs <- sort(nc$rhs)
        tmp <- paste(nc$season,prhs[1], sep=" -> ");
        if(length(prhs) > 1) for(i in 2:length(prhs)) tmp <- paste(tmp,prhs[i],sep=", ")
        return (tmp)
    }
            
    return (nc);
}

SeasonalFTS <- function(fsets,flrgs,pperiod){
    nc <- list(
        name = "Seasonal",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets),
        period = pperiod,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(f in seq(1,nc$period)){
            k <- nc$flrg[[f]]
            if(is.null(k)) k <- SFLRG(f, NULL)
            tmp <- sprintf("%s ; %s",tmp,k$dump());
        }
        return (tmp)
    }
    
    nc$getMidpoints <- function(fsets){
        mp <- c()
        
        if(length(fsets) == 0){
            return (0)
        }
        
        for(i in 1:length(fsets)) mp[i] <- nc$fuzzySets[[ fsets[i] ]]$midpoint;
        return (mp);
    }
    
    nc$forecast <- function(season){
		l <- length(season)
		ret <- c()

		for(k in 1:l) {
			fsets <- nc$flrg[[ season[k] ]];
			mp <- nc$getMidpoints( fsets$rhs )
			ret[k] <- (sum(mp)/length(mp))
			
		}	
        return ( ret )       
        
    }
    
    nc$forecastAhead <- function(x,steps){
		ret <- c()
		
		idx <- c()

		# Pay attention on the index correction of trailling values when mod = 0

		for(k in 1:(steps + (steps %/% (nc$period + 1)) )) 
			if(((x+k-1) %% (nc$period + 1)) > 0) 
				idx[k - ((x+k-1) %/% (nc$period + 1))] <- ((x+k-1) %% (nc$period + 1)) 
		
		ret <- nc$forecast(idx)
			
        return ( ret )       
    }
    
    return (nc)
}



FitSeasonalFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A"),
        period = parameters,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
    nc$genFLRG <- function(flrs){
        flrgs <- list()
        season <- 1
        for(flr in flrs){
            if(length(flrgs) < season){
                flrgs[[season]] <- SFLRG(season,c(flr$rhs));
            } else {
                flrgs[[season]] <- flrgs[[season]]$put(flr$rhs);
            }
                        
            season <- (season + 1) %% (nc$period + 1)
            if(season == 0) season <- 1
        }
        return (flrgs)
    } 
    
    nc$train <- function() {
       fzydata <- fuzzySeries(nc$data,nc$fuzzySets);
        flrs <- genRecurrentFLR(fzydata);
        flrgs <- nc$genFLRG(flrs);
        tmp <- SeasonalFTS(nc$fuzzySets, flrgs, nc$period);
        return (tmp);
    }
    
    return (nc)
}
