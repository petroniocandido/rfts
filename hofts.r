
HOFLRG <- function(plhs, prhs){
    nc <- list(
        lhs = plhs,
        rhs = prhs
    );
    
    nc$putLHS = function(x){
        if(!(x %in% nc$lhs)){
            return (HOFLRG(c(nc$lhs,x), nc$rhs));            
        } else { 
            return (nc);
        }
    }
    
    nc$putRHS = function(x){
        if(!(x %in% nc$rhs)){
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
            tmp <- sprintf("%s \n %s",tmp,k$dump());
        }
        return (tmp)
    }
    
    nc$getMidpoints <- function(nflrg){
        k <- nc$flrg[[ nflrg ]];
        mp <- c()
        
        if(length(k$rhs) == 0) {
			print("inexistente")
            return (c( nc$fuzzySets[[ nflrg ]]$midpoint))
        }
        
        for(i in 1:length(k$rhs)) mp[i] <- nc$fuzzySets[[ k$rhs[i] ]]$midpoint;
        return (mp);
    }
    
    nc$forecast <- function(x){
		l <- length(x)
		
		if(l < nc$lags){
			print("Dados insuficientes!")
		}

		ret <- c()

		for(k in nc$lags:l) {
			
			prev <- fuzzySeries( x[ (k - nc$lags):k ],  nc$fuzzySets)
			
			lhs <- toString(prev)
			
			mp <- nc$getMidpoints( lhs )
			
			ret[k] <- (sum(mp)/length(mp))
			
		}	
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
        for(k in nc$lags:l){
			
			prev <- fuzzySeries( data[ (k - nc$lags):k ],  nc$fuzzySets)
			
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
