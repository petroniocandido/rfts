SongFTS <- function(fsets,opm){
    nc <- list(
        name = "Song & Chissom",
        fuzzySets = fsets,
        operationMatrix = opm,
        npart = length(fsets)
    );
    
    nc$forecast <- function(x) {
        mv <- c();
        for(i in 1:nc$npart) mv[i] <- nc$fuzzySets[[i]]$membership(x);
        rmin <- matrix(rep(0,nc$npart*nc$npart), nrow=nc$npart,ncol=nc$npart);
        rmax <- c();
        for(i in 1:nc$npart) for(j in 1:nc$npart) rmin[i,j] <- min( nc$operationMatrix[i,j], mv[j] );
        for(i in 1:nc$npart) rmax[i] <- max( rmin[i,] );
        
        best_sets <- which(rmax == max(rmax) );
        
        if(length(best_sets) == 1){
            return (nc$fuzzySets[[best_sets[1]]]$midpoint);
        } else {
            vals <- c()
            for(i in 1:length(best_sets)) vals[i] <- nc$fuzzySets[[best_sets[i]]]$midpoint;
                           
            return ( sum(vals)/length(vals) )    
        } 
    }
    
            
    return (nc)
}

FitSongFTS <- function(pdata,np,mf,parameters) {
    nc <- list(
        data = pdata,
        npart = np,
        membershipFunc = mf,
        fuzzySets = UniversePartitioner(pdata,np,mf,"A")
    )
    
    nc$FLRmembershipMatrix <- function(pflr){
        ls <- nc$fuzzySets[[pflr$lhs]];
        lm <- c();
        rs <- nc$fuzzySets[[pflr$rhs]];
        rm <- c();
        mv <- matrix(nrow=nc$npart,ncol=nc$npart);
        for(i in 1:nc$npart){
            lm[i] <- ls$membership( nc$fuzzySets[[i]]$midpoint );
            rm[i] <- rs$membership( nc$fuzzySets[[i]]$midpoint );
        }
        
        for(i in 1:nc$npart) for(j in 1:nc$npart) mv[i,j] <- min(lm[j],rm[i]);
        
        return (mv);
    }
    
    nc$operationMatrix <- function(flrs) {
        r  <- matrix( rep(0,nc$npart*nc$npart), nrow=nc$npart,ncol=nc$npart);
        for(k in 1:length(flrs)){
            mm <- nc$FLRmembershipMatrix( flrs[[k]] );
            for(i in 1:nc$npart) for(j in 1:nc$npart) r[i,j] <- max( r[i,j], mm[i,j] );
        }
        return (r);
    }
    
    nc$train <- function() {
        fzydata <- fuzzySeries(nc$data,nc$fuzzySets);
        flrs <- genFLR(fzydata);
        r <- nc$operationMatrix(flrs);
        tmp <- SongFTS(nc$fuzzySets, r);
        return (tmp);
    }
    return (nc);
}
