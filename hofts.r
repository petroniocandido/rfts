
HOFLRG <- function(plhs, prhs){
    nc <- list(
        lhs = plhs,
        rhs = prhs
    );
    
    nc$putLHS = function(x){
        if(!(x %in% nc$lhs) && !is.na(x) ){
            return (HOFLRG(c(nc$lhs,x), nc$rhs));            
        } else { 
            return (nc);
        }
    }
    
    nc$putRHS = function(x){
        if(!(x %in% nc$rhs) && !is.na(x) ){
            return (HOFLRG(nc$lhs, c(nc$rhs,x)));            
        } else { 
            return (nc);
        }
    } 
            
    nc$dump <- function() {
		tmpl <- toString(nc$lhs)
		tmpr <- toString(sort(nc$rhs)) 
        tmp <- paste(tmpl,tmpr, sep=" -> ");
        return (tmp)
    }
            
    return (nc);
}

HOFTS <- function(fsets,flrgs,l){
    nc <- list(
        name = "HOFTS",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets),
        lags = l,
        isHighOrder = TRUE,
        isIntervallic = FALSE
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(k in nc$flrg){
            if(is.null(k)) k <- HOFLRG(fs$name, c(fs$name))
            tmp <- sprintf("%s ; %s",tmp,k$dump());
        }
        return (tmp)
    }
    
    nc$getMidpoints <- function(nflrg){
        k <- nc$flrg[[ nflrg ]];
        mp <- c()
        if(length(k$rhs) == 0) {
			tmp <- unlist(strsplit(nflrg,", "))
			ltmp <- length(tmp)
			mp <- c(nc$fuzzySets[[ tmp[ltmp] ]]$midpoint);
			
			return (mp);
        } else {
			for(i in 1:length(k$rhs)) mp[i] <- nc$fuzzySets[[ k$rhs[i] ]]$midpoint;
			return (mp);
		}
    }
    
    nc$forecast <- function(x){
		
		x <- x[ !is.na(x) ]
		l <- length(x)
		
		if(l < nc$lags){
			return (c(x))
		}

		ret <- c()
		c <- 1
		
		for(k in seq(nc$lags,l)) {
			prev <- fuzzySeries( x[ seq(k + 1 - nc$lags,k) ],  nc$fuzzySets)
			lhs <- toString(prev)
			mp <- nc$getMidpoints( lhs )
			ret[c] <- as.numeric(sum(mp)/as.numeric(length(mp)))
			c <- c + 1
		}	
				
        return ( ret )       
        
    }
    
    nc$forecastAhead <- function(x,steps){
		ret <- x

		for(k in nc$lags:steps) ret[k+1] <- nc$forecast(ret[seq(k-nc$lags,k)])
		
        return ( ret )       
    }
    
    return (nc)
}


FitHOFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A"),
        lags = parameters,
        isHighOrder = TRUE,
        isIntervallic = FALSE
    );
    
    nc$genFLRG <- function(data){
        flrgs <- list()
        l <- length(data)
        for(k in nc$lags+1:l){
			
			prev <- data[ seq(k - nc$lags,k-1) ]
			
			lhs <- toString(prev)
			
			if(is.null(flrgs[[lhs]])){
                flrgs[[lhs]] <- HOFLRG(prev,c(data[k]));
            } else {
                flrgs[[lhs]] <- flrgs[[lhs]]$putRHS(data[k]);
            }
        }
        return (flrgs)
    } 
    
    nc$train <- function() {
       fzydata <- fuzzySeries(nc$data,nc$fuzzySets);
        flrgs <- nc$genFLRG(fzydata);
        tmp <- HOFTS(nc$fuzzySets, flrgs, nc$lags);
        return (tmp);
    }
    
    return (nc)
}
