D <- function(data, t){
    newdata <- c()
    for(i in 1:length(data)-t){
        newdata[i] <- data[i+1] - data[i]
    }
    return (newdata)
}

trimf <- function(x,parameters){
    if(x < parameters[1]){
		#print("antes")
        return (0)
    }
    else if(x >= parameters[1] && x < parameters[2]){
		#print("subindo")
        return ((x-parameters[1])/(parameters[2]-parameters[1]))
    }
    else if(x >= parameters[2] && x <= parameters[3]) {
		#print("descendo")
        return ((parameters[3]-x)/(parameters[3]-parameters[2]))
    }
    else { 
		#print("depois")
        return (0)
    }
}
    
FuzzySet <- function(pname, pmf, pparameters,pmidpoint){
    nc = list(
        name = pname,
        parameters = pparameters,
        lower = min(pparameters),
        upper = max(pparameters),
        midpoint = pmidpoint,
        mf = pmf
    )
    nc$membership = function(x){ 	
		return (nc$mf(x,nc$parameters)) 
    }
    nc$dump = function(){ 
        pars <-paste(sprintf("%s", nc$parameters) )
        print(sprintf("%s : %s",nc$name, pars)) 
    }
    return (nc)
}
