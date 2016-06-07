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
        return (min(lw));
    }
        
    nc$getUpper <- function(nflrg){
        k <- nc$flrg[[nflrg]];
         if(is.null(k))
            return (nc$fuzzySets[[ nflrg ]]$upper)
             
        if(length(k$rhs) == 0)
            return (nc$fuzzySets[[ nflrg ]]$upper)
            
        lw <- c()
        for(i in 1:length(k$rhs)) lw[i] <- nc$fuzzySets[[ k$rhs[i] ]]$upper;
        return (max(lw));
    }
    
    nc$forecast <- function(x){
        mv <- c(); 
        lw <- c();
        up <- c();
        #print(x)
        for(i in 1:nc$npart) { 
            fs <- nc$fuzzySets[[i]];
            mv[i] <- fs$membership(x);
            #print(fs$membership)
            lw[i] <- mv[i] * nc$getLower( fs$name );
            up[i] <- mv[i] * nc$getUpper( fs$name );
        }
        
        return ( matrix(c( sum(lw), sum(up) ), 1,2) )
        
    }
        
    nc$forecastAhead <- function(x,n){
        forecasts <- matrix( rep(0,2*n), n,2)
        interval <- matrix( rep(0,20), 10,2)
        
        forecasts[1,] <- nc$forecast(x)
        for(i in 2:n){
            inc <- (forecasts[i-1,2] - forecasts[i-1,1])/9
            c <- 1
            for(k in seq(forecasts[i-1,1], forecasts[i-1,2], inc)){
                forecasts[i,] <- forecasts[i,] + nc$forecast(k)
                c <- c + 1
            }
            forecasts[i,] = forecasts[i,] / (c-1)
        }
        return (forecasts)
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
