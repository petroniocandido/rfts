benchmark <- function(builder,pdata,np,mf,parameters,crossvalidation_ratio){
    size <- length(pdata)
    train <- pdata[1:round(size*crossvalidation_ratio)]
    test <- pdata[round(size*crossvalidation_ratio):size]

    tmp <- builder(train,np,mf,parameters)
    
    model <- tmp$train()
    
    tryCatch(pred <- sapply(test, model$forecast), error= function(e){ print(e); print(sprintf(model$dump())) })
    
    return (list(model = model, 
                 validation = c(test,NA), predicted = c(NA,pred), 
                 mape = MAPE(c(test,NA),c(NA,pred)),
                 rmse = RMSE(c(test,NA),c(NA,pred)),
                 validationIndex = round(size*crossvalidation_ratio)))
}

benchmarkInterval <- function(builder,pdata,np,mf,parameters,crossvalidation_ratio){
    size <- length(pdata)
    train <- pdata[1:round(size*crossvalidation_ratio)]
    test <- pdata[round(size*crossvalidation_ratio):size]

    tmp <- builder(train,np,mf,parameters)
    
    model <- tmp$train()
    
    tryCatch(pred <- sapply(test, model$forecast), error= function(e){ print(e); print(sprintf(model$dump())) })
    
    return (list(model = model, 
                 validation = c(test,NA), predicted = c(NA,pred), 
                 mape = MAPE(c(test,NA),c(NA,pred)),
                 rmse = RMSE(c(test,NA),c(NA,pred)),
                 validationIndex = round(size*crossvalidation_ratio)))
}

#
# Train a model in many partitions sizes, checking the error mesures
#

testPartitions <- function(builder,pdata,indexfield,valuefield,nps,crossvalidation_ratio){
    predictions <- data.frame(validation<-c(),predicted<-c(),model<-c())
    rmse <- c()
    mape <- c()
    benchmarks <- list()
    for(i in 1:length(nps)){
        benchmarks[[i]] <- benchmark(builder, as.vector(pdata[,valuefield]),nps[i],
                                   trimf,NULL,crossvalidation_ratio)
        tmp <- data.frame( benchmarks[[i]]$validation,
                          benchmarks[[i]]$predicted,
                          rep(nps[i],length(benchmarks[[i]]$predicted)) )
        names(tmp) <- c("validation","predicted","model")
        predictions <- rbind(predictions, tmp)
        rmse[i] <- benchmarks[[i]]$rmse
        mape[i] <- benchmarks[[i]]$mape
    }
    options(repr.plot.width=9, repr.plot.height=7)
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    boxplot(abs(validation-predicted) ~ model, data=predictions, main="Error distribution by number of partitions",
           xlab="Partitions",ylab="Absolute error")
    plot(nps,rmse,xlab="Partitions",ylab="Error",main="RMSE",type="o")
    plot(nps,mape,xlab="Partitions",ylab="Error",main="MAPE",type="o")
}


#
# Perform all validation tasks in all models for a dataset
#

benchmarkAll <- function(data,indexfield,valuefield,nps,crossvalidation_ratio){
    builders <- list(Song = FitSongFTS, Chen = FitChenFTS, Yu = FitYuFTS)
    benchmarks <- list()
    for(i in 1:length(builders)) benchmarks[[i]] <- benchmark(builders[[i]],
                                                              as.vector(data[,valuefield]),nps[i],trimf,NULL,
                                                              crossvalidation_ratio)
    td <- length(data[,indexfield])
    vh <- mean(data[,valuefield])
    vi <- benchmarks[[1]]$validationIndex
    
    options(repr.plot.width=9, repr.plot.height=10)
    layout(matrix(c(1,1,1,2,2,2,3,4,5), 3, 3, byrow = TRUE))
        
    plot(data[,indexfield],data[,valuefield],
         type="l",main="Original data and cross validation sets",
        xlab="t",ylab="F(t)")
    text(min(data[,indexfield]),vh,"Training set",pos=4)
    text(data[vi,indexfield],vh,"Validation set",pos=4)
    abline(h=vh,v=data[vi,indexfield])
    
    test <- benchmarks[[1]]$validation
        
    miny <- min(test,na.rm=TRUE)
    maxy <- max(test,na.rm=TRUE)    
    for(i in 1:length(builders)) {
        miny <- min(miny, min(benchmarks[[i]]$predicted,na.rm=TRUE))
        maxy <- max(maxy, max(benchmarks[[i]]$predicted,na.rm=TRUE))
    }
        
    plot(test, type="l", col=1, 
         main="Model performance on validation set",xlab="t",ylab="F(t)",
         ylim=c(miny*0.9,maxy*1.1))
    lgd <- c("Original")
    clrs <- c(1)
    for(i in 1:length(builders)) {
        lgd[i+1] <- benchmarks[[i]]$model$name 
        clrs[i+1] <- i*7
        lines(benchmarks[[i]]$predicted,type="l",col=i*7)
    }
    legend("topright",legend=lgd, fill=clrs)
    
    predictions <- data.frame(validation<-c(),predicted<-c(),model<-c())
    rmse <- c()
    mape <- c()
    mdl <- c()
    for(i in 1:length(builders)) {
        tmp <- data.frame( benchmarks[[i]]$validation,
                          benchmarks[[i]]$predicted,
                          rep(benchmarks[[i]]$model$name,length(benchmarks[[i]]$predicted)) )
        names(tmp) <- c("validation","predicted","model")
        predictions <- rbind(predictions, tmp)
        rmse[i] <- benchmarks[[i]]$rmse
        mape[i] <- benchmarks[[i]]$mape
        mdl[i] <- benchmarks[[i]]$model$name
    }
    boxplot(abs(validation-predicted) ~ model, data=predictions, main="Error distribution by model")
    
    tmp <- data.frame(rmse,mape,mdl)
        
    plot(rmse ~ mdl, data=tmp, main="RMSE",type="l",lwd=5,ylab="",xlab="")
    plot(mape ~ mdl, data=tmp, main="MAPE",type="h",lwd=5,ylab="",xlab="")
}


