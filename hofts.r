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
		tmpl <- lhs[1]
		if(length(lhs) > 1) for(i in 2:length(lhs)) tmpl <- paste(tmpl,lhs[i],sep=", ")
        prhs <- sort(nc$rhs)
        tmp <- paste(tmpl,prhs[1], sep=" -> ");
        if(length(prhs) > 1) for(i in 2:length(prhs)) tmp <- paste(tmp,prhs[i],sep=", ")
        return (tmp)
    }
            
    return (nc);
}
