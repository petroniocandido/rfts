EWFLRG <- function(plhs, prhs, pc){
    nc <- list(
        lhs = plhs,
        rhs = prhs,
        c = as.numeric(pc)
    );
    
    nc$put = function(x){
       return (EWFLRG(nc$lhs, c(nc$rhs,x),nc$c));       
    }
    
    nc$getTotalWeight <- function(){
		tot <- 0.0
		for(i in 1:length(nc$rhs)){
			tot <- tot + as.numeric(nc$c**(i-1))
		}
		return (tot)
    }
    
    nc$getWeight <- function(x){
		tot <- nc$getTotalWeight()
		return ((nc$c**(x-1))/tot)
    }
    
    nc$getWeights = function(){
		
		 w <- c();
        
        if(length(nc$rhs) == 0)
            return (matrix(c(1)))
                            
        for(i in 1:length(nc$rhs)) w[i] <- as.numeric(nc$c**(i-1));
        
        tot <- sum(w)
        
        w <- w / tot
        
        #print(w)
        
        return (matrix(w))
		
	}
	
    nc$dump <- function() {
        prhs <- nc$rhs
        tmp <- paste(nc$lhs, sprintf("%s*%s",round(nc$getWeight(1),2),nc$rhs[1]), sep=" -> ");
        
        if(length(prhs) > 1) 
            for(i in 2:length(nc$rhs)) 
                tmp <- paste(tmp,sprintf("%s*%s",round(nc$getWeight(i),2),nc$rhs[i]),sep=", ")
        return (tmp)
    }
            
    return (nc);
}

SadaeiFTS <- function(fsets,flrgs,pc){
    nc <- list(
        name = "Sadaei",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets),
        c = pc,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(fs in nc$fuzzySets){
            k <- nc$flrg[[ fs$name ]]
            if(is.null(k)) k <- EWFLRG(fs$name, c(fs$name),nc$c)
            tmp <- sprintf("%s ; %s",tmp,k$dump());
        }
        return (tmp)
    }
    
    nc$getMidpoints <- function(nflrg){
        
        k <- nc$flrg[[nflrg]];
        mp <- c()
        
        if(is.null(k))
            return (matrix(c(nc$fuzzySets[[ nflrg ]]$midpoint)))
            
        if(length(k$rhs) == 0)
            return (matrix(c(nc$fuzzySets[[ nflrg ]]$midpoint)))
        
        for(i in 1:length(k$rhs)) mp[i] <- nc$fuzzySets[[ k$rhs[i] ]]$midpoint;
        return (matrix(mp));
    }
    
    nc$forecast <- function(x){
		
		l <- length(x)

		ret <- c()

		for(k in 1:l) {
			mv <- c();
        
			for(i in 1:nc$npart) mv[i] <- nc$fuzzySets[[i]]$membership(x[k]);
			
			best_sets <- which(mv == max(mv));
			
			lhs <- nc$fuzzySets[[ best_sets[1] ]]$name
			
			mp <- nc$getMidpoints( lhs )
						
			wg <- matrix(c(1))
			if(lhs %in% names(nc$flrg)){
				wg <- nc$flrg[[lhs]]$getWeights()
			}
			
			
			ret[k] <- (t(wg) %*% mp)
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

FitSadaeiFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A"),
        c = parameters,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
    nc$genFLRG <- function(flrs){
        flrgs <- list()
        for(flr in flrs){
            #print(flr$dump())
            if(is.null(flrgs[[flr$lhs]])){
                flrgs[[flr$lhs]] <- EWFLRG(flr$lhs,c(flr$rhs),nc$c);
            } else {
                flrgs[[flr$lhs]] <- flrgs[[flr$lhs]]$put(flr$rhs);
            }
        }
        return (flrgs)
    } 
    
    nc$train <- function() {
       fzydata <- fuzzySeries(nc$data,nc$fuzzySets);
       flrs <- genRecurrentFLR(fzydata);
       flrgs <- nc$genFLRG(flrs);
       tmp <- SadaeiFTS(nc$fuzzySets, flrgs,nc$c);
       return (tmp);
    }
    
    return (nc)
}
