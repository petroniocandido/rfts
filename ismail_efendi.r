IWFLRG <- function(plhs, prhs, ptotal){
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
       return (IWFLRG(nc$lhs, nc$rhs, nc$total));       
    }
    
    nc$getWeight = function(x){
		if(is.null(nc$rhs[[x]])){
			return (0)
		} else {
			return(nc$rhs[[x]]/nc$total)
		}
	}
	
	nc$getWeights = function(){
		
		 w <- c();
        
        if(length(nc$rhs) == 0)
            return (matrix(c(1)))
        
        for(i in 1:length(nc$rhs)) w[i] <- nc$rhs[[i]]/nc$total;
        return (matrix(w))
		
	}
	
    nc$dump <- function() {
        prhs <- nc$rhs
        tmp <- paste(sprintf("\n"), paste(nc$lhs, sprintf("%s*%s",round(nc$getWeight(names(nc$rhs)[1]),2),names(nc$rhs)[1]), sep=" -> "));
        if(length(prhs) > 1) 
            for(i in 2:length(nc$rhs)) 
                tmp <- paste(tmp,sprintf("%s*%s",round(nc$getWeight(names(nc$rhs)[i]),2),names(nc$rhs)[i]),sep=", ")
        return (tmp)
    }
            
    return (nc);
}

EfendiFTS <- function(fsets,flrgs){
    nc <- list(
        name = "Ismail & Efendi",
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
            if(is.null(k)) {
				 tmp <- list()
				 tmp[[fs$name]] <- 1
				 k <- IWFLRG(fs$name, tmp, 1)
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
        
        for(i in 1:length(k$rhs)) {
			mp[i] <- nc$fuzzySets[[ names(k$rhs)[i] ]]$midpoint;
		}
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

FitEfendiFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A"),
        name = "Improve Weighted FTS",
        isSeasonal = FALSE,
        isHighOrder = FALSE,
        isIntervallic = FALSE
    );
    
     nc$genFLRG <- function(flrs){
        flrgs <- list()
        for(flr in flrs){
            if(is.null(flrgs[[flr$lhs]])){
				tmp <- list()
				tmp[[flr$rhs]] = 1
                flrgs[[flr$lhs]] <- IWFLRG(flr$lhs,tmp,1);
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
       tmp <- EfendiFTS(nc$fuzzySets, flrgs);
       return (tmp);
    }
    
    return (nc)
}
