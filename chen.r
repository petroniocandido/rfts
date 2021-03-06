
FLRG <- function(plhs, prhs){
    nc <- list(
        lhs = plhs,
        rhs = prhs
    );
    
    nc$put = function(x){
        if(!(x %in% nc$rhs)){
            return (FLRG(nc$lhs, c(nc$rhs,x)));            
        } else { 
            return (nc);
        }
    } 
            
    nc$dump <- function() {
        prhs <- sort(nc$rhs)
        tmp <- paste(nc$lhs,prhs[1], sep=" -> ");
        if(length(prhs) > 1) for(i in 2:length(prhs)) tmp <- paste(tmp,prhs[i],sep=", ")
        return (tmp)
    }
            
    return (nc);
}

ChenFTS <- function(fsets,flrgs){
    nc <- list(
        name = "Chen",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets),
        isSeasonal = FALSE,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(fs in nc$fuzzySets){
            k <- nc$flrg[[ fs$name ]]
            if(is.null(k)) k <- FLRG(fs$name, c(fs$name))
            tmp <- sprintf("%s \n %s",tmp,k$dump());
        }
        return (tmp)
    }
    
    nc$getMidpoints <- function(nflrg){
        k <- nc$flrg[[ nflrg ]];
        mp <- c()
        
        if(length(k$rhs) == 0)
            return (c( nc$fuzzySets[[ nflrg ]]$midpoint))
        
        for(i in 1:length(k$rhs)) mp[i] <- nc$fuzzySets[[ k$rhs[i] ]]$midpoint;
        return (mp);
    }
    
    nc$forecast <- function(x){
		l <- length(x)

		ret <- c()

		for(k in 1:l) {
			mv <- c();
			for(i in 1:nc$npart) mv[i] <- nc$fuzzySets[[i]]$membership(x[k]);
			
			best_sets <- which(mv == max(mv));
			
			mp <- nc$getMidpoints( nc$fuzzySets[[ best_sets[1] ]]$name )
			
			ret[k] <- (sum(mp)/length(mp))
			
		}	
        return ( ret )       
        
    }
    
    nc$forecastAhead <- function(x,steps){
		ret <- c(x)

		for(k in 1:steps) {
			ret[k+1] <- nc$forecast(ret[k])
		}	
        return ( ret )       
    }
    
    return (nc)
}



FitChenFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A"),
        name = "Conventional FTS",
        isSeasonal = FALSE,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
    nc$genFLRG <- function(flrs){
        flrgs <- list()
        for(flr in flrs){
            if(is.null(flrgs[[flr$lhs]])){
                flrgs[[flr$lhs]] <- FLRG(flr$lhs,c(flr$rhs));
            } else {
                flrgs[[flr$lhs]] <- flrgs[[flr$lhs]]$put(flr$rhs);
            }
        }
        return (flrgs)
    } 
    
    nc$train <- function() {
       fzydata <- fuzzySeries(nc$data,nc$fuzzySets);
        flrs <- genFLR(fzydata);
        flrgs <- nc$genFLRG(flrs);
        tmp <- ChenFTS(nc$fuzzySets, flrgs);
        return (tmp);
    }
    
    return (nc)
}
