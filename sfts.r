SLRG <- function(plhs, prhs){
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
