PFLRG <- function(plhs, prhs){
    nc <- list(
        lhs = plhs,
        rhs = prhs
    );
    
    nc$put = function(x){
       ind = which(nc$rhs == x)
       return (PFLRG(nc$lhs, c(nc$rhs,x)));       
    }
    
    nc$dump <- function() {
        prhs <- sort(nc$rhs)
        tmp <- paste(nc$lhs, prhs[1], sep=" -> ");
        if(length(prhs) > 1) 
            for(i in 2:length(prhs)) 
                tmp <- paste(tmp,prhs[i],sep=", ")
        return (tmp)
    }
            
    return (nc);
}

PWFLRG <- function(plhs, prhs, ptotal){
    nc <- list(
        lhs = plhs,
        rhs = prhs,
        total = ptotal
    );
    
    nc$put = function(x){
       if(x %in% names(nc$rhs)){
		   nc$rhs[[x]] <- nc$rhs[[x]] + 1
	   } else {
		   nc$rhs[[x]] <- 1
	   }
	   nc$total <- nc$total + 1
       return (PWFLRG(nc$lhs, nc$rhs, nc$total));       
    }
    
    nc$getWeight = function(x){
		if(is.null(nc$rhs[[x]])){
			return (0)
		} else {
			return(nc$rhs[[x]]/nc$total)
		}
	}
	
    nc$dump <- function() {
        prhs <- nc$rhs
        tmp <- paste(sprintf("\n"), paste(nc$lhs, sprintf("%s*%s",nc$getWeight(names(nc$rhs)[1]),names(nc$rhs)[1]), sep=" -> "));
        if(length(prhs) > 1) 
            for(i in 2:length(nc$rhs)) 
                tmp <- paste(tmp,sprintf("%s*%s",nc$getWeight(names(nc$rhs)[i]),names(nc$rhs)[i]),sep=", ")
        return (tmp)
    }
            
    return (nc);
}

PFTS <- function(fsets,flrgs){
    nc <- list(
        name = "Probabilistic FTS",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets)
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(fs in nc$fuzzySets){
            k <- nc$flrg[[ fs$name ]]
            if(is.null(k)) k <- PFLRG(fs$name, c(fs$name))
            tmp <- sprintf("%s \n %s",tmp,k$dump());
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
        
    nc$getLower <- function(nflrg){
        k <- nc$flrg[[nflrg]];
         if(is.null(k))
            return (nc$fuzzySets[[ nflrg ]]$lower)
             
        if(length(k$rhs) == 0)
            return (nc$fuzzySets[[ nflrg ]]$lower)
            
        lw <- c()
        for(i in 1:length(k$rhs)) lw[i] <- nc$fuzzySets[[ k$rhs[i] ]]$lower;
        return (mean(lw));
    }
        
    nc$getUpper <- function(nflrg){
        k <- nc$flrg[[nflrg]];
         if(is.null(k))
            return (nc$fuzzySets[[ nflrg ]]$upper)
             
        if(length(k$rhs) == 0)
            return (nc$fuzzySets[[ nflrg ]]$upper)
            
        lw <- c()
        for(i in 1:length(k$rhs)) lw[i] <- nc$fuzzySets[[ k$rhs[i] ]]$upper;
        return (mean(lw));
    }
    
    nc$forecast <- function(x){
		l <- length(x)

		ret <- matrix(rep(0,l*2), l,2)

		for(k in 1:l) {
			mv <- c(); 
			lw <- c();
			up <- c();
			for(i in 1:nc$npart) { 
				fs <- nc$fuzzySets[[i]];
				mv[i] <- fs$membership(x[k]);
				lw[i] <- mv[i] * nc$getLower( fs$name );
				up[i] <- mv[i] * nc$getUpper( fs$name );
			}
			ret[k,] <- c( sum(lw), sum(up) )
		}
        return ( ret )
    }
        
    nc$forecastAhead <- function(x,n){
        forecasts <- matrix( rep(0,2*n), n,2)
        interval <- matrix( rep(0,20), 10,2)
        
        forecasts[1,] <- nc$forecast(x)
        for(i in 2:n){
			lower <- nc$forecast(forecasts[i-1,1])
			upper <- nc$forecast(forecasts[i-1,2])
			forecasts[i,1] <- lower[1,1]
			forecasts[i,2] <- upper[1,2]
            #inc <- (forecasts[i-1,2] - forecasts[i-1,1])/10
            #c <- 1
            #for(k in seq(forecasts[i-1,1], forecasts[i-1,2], inc)){
				#tmp <- nc$forecast(k)
                #forecasts[i,1] <- forecasts[i,1] + tmp[1,1]
                #forecasts[i,1] <- min(forecasts[i,1],tmp[1,1])
                #forecasts[i,2] <- forecasts[i,2] + tmp[1,2]
                #forecasts[i,2] <- max(forecasts[i,2],tmp[1,2])
                #c <- c + 1
            #}
            #forecasts[i,] = forecasts[i,] / (c-1)
        }
        return (forecasts)
    }
    
    return (nc)
}

PWFTS <- function(fsets,flrgs){
    nc <- list(
        name = "Probabilistic Weighted FTS",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets)
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(fs in nc$fuzzySets){
            k <- nc$flrg[[ fs$name ]]
            if(is.null(k)) {
				 tmp <- list()
				 tmp[[fs$name]] <- 1
				 k <- PWFLRG(fs$name, tmp, 1)
			}
            tmp <- sprintf("%s \n %s",tmp,k$dump());
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
        
        for(i in 1:length(k$rhs)) mp[i] <- nc$fuzzySets[[ k$rhs[[i]] ]]$midpoint;
        return (matrix(mp));
    }
        
    nc$getLower <- function(nflrg){
        k <- nc$flrg[[nflrg]];
         if(is.null(k))
            return (nc$fuzzySets[[ nflrg ]]$lower)
             
        if(length(k$rhs) == 0)
            return (nc$fuzzySets[[ nflrg ]]$lower)
            
        lw <- 0
        for(i in names(k$rhs)) lw <- lw + k$getWeight(i) * nc$fuzzySets[[ i ]]$lower;
        return (lw);
    }
        
    nc$getUpper <- function(nflrg){
		
        k <- nc$flrg[[nflrg]];
         if(is.null(k))
            return (nc$fuzzySets[[ nflrg ]]$upper)
             
        if(length(k$rhs) == 0)
            return (nc$fuzzySets[[ nflrg ]]$upper)
            
        lw <- 0
        for(i in names(k$rhs)) lw <- lw + k$getWeight(i) * nc$fuzzySets[[ i ]]$upper;
        
        return (lw);
    }
    
    nc$forecast <- function(x){
		l <- length(x)

		ret <- matrix(rep(0,l*2), l,2)

		for(k in 1:l) {
			mv <- c(); 
			lw <- c();
			up <- c();
			for(i in 1:nc$npart) { 
				fs <- nc$fuzzySets[[i]];
				mv[i] <- fs$membership(x[k]);
				lw[i] <- mv[i] * nc$getLower( fs$name );
				up[i] <- mv[i] * nc$getUpper( fs$name );
			}
			ret[k,] <- c( sum(lw), sum(up) )
		}
        return ( ret )
    }
    
    
    return (nc)
}


FitPFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A")
    );
    
    nc$genFLRG <- function(flrs){
        flrgs <- list()
        for(flr in flrs){
            if(is.null(flrgs[[flr$lhs]])){
                flrgs[[flr$lhs]] <- PFLRG(flr$lhs,c(flr$rhs));
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
       tmp <- PFTS(nc$fuzzySets, flrgs);
       return (tmp);
    }
    
    return (nc)
}


FitPWFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A")
    );
    
    nc$genFLRG <- function(flrs){
        flrgs <- list()
        for(flr in flrs){
            if(is.null(flrgs[[flr$lhs]])){
				tmp <- list()
				tmp[[flr$rhs]] = 1
                flrgs[[flr$lhs]] <- PWFLRG(flr$lhs,tmp,1);
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
       tmp <- PWFTS(nc$fuzzySets, flrgs);
       return (tmp);
    }
    
    return (nc)
}
