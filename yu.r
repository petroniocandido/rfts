YuFLRG <- function(plhs, prhs){
    nc <- list(
        lhs = plhs,
        rhs = prhs
    );
    
    nc$put = function(x){
       return (YuFLRG(nc$lhs, c(nc$rhs,x)));       
    }
    
    nc$getWeight <- function(x){
        return (x/(sum(seq(1,length(nc$rhs)))))
    }
    
    nc$getWeights <- function(){
        w <- c();
        
        if(length(nc$rhs) == 0)
            return (matrix(c(1)))
        
        len <- sum(seq(1,length(nc$rhs)));
        for(i in 1:length(nc$rhs)) w[i] <- i/len;
        return (matrix(w))            
    }
            
    nc$dump <- function() {
        prhs <- sort(nc$rhs)
        tmp <- paste(nc$lhs, sprintf("%s * %s",round(nc$getWeight(1),2),prhs[1]), sep=" -> ");
        if(length(prhs) > 1) 
            for(i in 2:length(prhs)) 
                tmp <- paste(tmp,sprintf("%s * %s",round(nc$getWeight(i),2),prhs[i]),sep=", ")
        return (tmp)
    }
            
    return (nc);
}

YuFTS <- function(fsets,flrgs){
    nc <- list(
        name = "Yu",
        fuzzySets = fsets,
        flrg = flrgs,
        npart = length(fsets)
    );
    
    nc$dump <- function() {
        tmp <- ""
        for(fs in nc$fuzzySets){
            k <- nc$flrg[[ fs$name ]]
            if(is.null(k)) k <- YuFLRG(fs$name, c(fs$name))
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
    
    nc$forecast <- function(x){
        mv <- c();
        
        for(i in 1:nc$npart) mv[i] <- nc$fuzzySets[[i]]$membership(x);
        
        best_sets <- which(mv == max(mv));
        
        lhs <- nc$fuzzySets[[ best_sets[1] ]]$name
        
        mp <- nc$getMidpoints( lhs )
        
        wg <- matrix(c(1))
        if(lhs %in% names(nc$flrg)){
            wg <- nc$flrg[[lhs]]$getWeights()
        }

        return (t(wg) %*% mp)
        
    }
    
    return (nc)
}

FitYuFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A")
    );
    
    nc$genFLRG <- function(flrs){
        flrgs <- list()
        for(flr in flrs){
            #print(flr$dump())
            if(is.null(flrgs[[flr$lhs]])){
                flrgs[[flr$lhs]] <- YuFLRG(flr$lhs,c(flr$rhs));
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
       tmp <- YuFTS(nc$fuzzySets, flrgs);
       return (tmp);
    }
    
    return (nc)
}
