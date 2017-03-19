md.dist <- function(d1, d2){
    d = d1 - d2
    mean(sqrt(diag(d%*%t(d))))
}

mddl <-function(obs, pred){
    TSDecomposition::mddl(obs, pred, plot=FALSE)
}

mda <- function(obs, pred, m, d){
    at1 = embedd(obs,  m=m, d=d)
    at2 = embedd(pred, m=m, d=d)

    n.elem = nrow(at1)
    delta = nrow(at2) - nrow(at1)

    if(delta == 0){
        return(md.dist(at1, at2))
    } else {
        return(min(lapply(0:delta, function(x) md.dist(at1, at2[(1+x):(n.elem+x),]))))
    }
}
        
loadSeriesFile <- function(seriesFolder){
    #select all series file in the data folder
    seriesFile = list.files(path = seriesFolder, full.names = TRUE)
    seriesList = list()
    for(i in 1:length(seriesFile)){
        load(seriesFile[i])
        seriesList[[i]] = seriesObj
    }
    seriesList
}

exec <- function(series, F, param){
    result = c()
    tryCatch({
        result = F(series, param)
    }, warning = function(w){
    }, error = function(e){
    }, finally = {
    })
    if(is.null(result)) result = rep(0, length(series))
    result
}
        
foreachParam <- function(series, F, hyperparameters){
    t(apply(hyperparameters, 1, function(x) exec(series, F, x)))
}
          
evaluateResult <- function(obs, resultSeries, m, d){
    resultTable = matrix(ncol=3,nrow=nrow(resultSeries))
    resultTable[,1:2] = apply(resultSeries, 
                              1, 
                              function(pred) c(TSDecomposition::mddl(obs, pred, plot=FALSE), mda(obs, pred, m, d)))
    standMddl = (resultTable[,1] - mean(resultTable[,1]))/(sd(resultTable[,1]))
    standMda  = (resultTable[,2] - mean(resultTable[,2]))/(sd(resultTable[,2]))
    resultTable[,3] = sqrt(standMddl^2 + standMda^2)
    colnames(resultTable) = c('MDDL', 'MDA', 'Dist')
    resultTable
}
    
gridSearch <- function(F, hyperparameters, seriesList, modelFolder, techName, cores = 1, verbose = FALSE){

    registerDoMC(cores)

    resultTable = foreach(i=1:length(seriesList), .combine='rbind') %dopar% {
        #load series object from RData file
        seriesObj = seriesList[[i]]
        m = unlist(seriesObj$det.embDim)
        d = unlist(seriesObj$det.sepDim)
        
        #Cross-validation out-of-time
        l = ceiling(0.7*seriesObj$size)
        trainIdx = 1:(l-1)
        testIdx  = l:seriesObj$size
    
        #Grid Search for better hyperparameters
        if(verbose) cat(paste('ts:',i,'- hyperparameters tunning: begin\n'))
        startTime = Sys.time()
        resultSeries   = foreachParam(seriesObj$series[trainIdx], F, hyperparameters)
        endTime = Sys.time()
        if(verbose) cat(paste('ts:',i,'- hyperparameters tunning: done - ', endTime - startTime,'\n'))
        
        if(verbose) cat(paste('ts:',i,'- evaluation: begin','\n'))
        startTime = Sys.time()
        validTestIdx   = which(rowSums(resultSeries) > 0)
        resultTableAux = evaluateResult(seriesObj$det.series[trainIdx], resultSeries[validTestIdx,], m, d)
        bestParamIdx   = which.min(resultTableAux[,3])
        resultTableAux = cbind( rep(techName,nrow(resultTableAux)),
                                apply(hyperparameters[validTestIdx], 1, function(x) paste(x, collapse=",")),
                                resultTableAux
                              )
        endTime = Sys.time()
        if(verbose) cat(paste('ts:',i,'- evaluation: done - ', endTime - startTime,'\n'))
                          
        if(verbose) cat(paste('ts:',i,'- test series evaluation: begin','\n'))
        startTime = Sys.time()
        #run best model with test data and evaluate it
        bestP = hyperparameters[validTestIdx[bestParamIdx],]
        pred  = F(seriesObj$series[testIdx], bestP)
        obs   = seriesObj$det.series[testIdx]
        mddlR = mddl(obs, pred)
        mdaR  = mda (obs, pred, seriesObj$det.embDim, seriesObj$det.sepDim)
        
        #save model result into model folder
        model = list( model.name = tech.name, F = F, best.param = bestP, mda = mdaR, mddl = mddlR )
        save(model, file=paste(modelFolder, '/', techName, '_', gsub("[^\\d]+", "", i, perl=TRUE) ,'.RData',sep=''))
        endTime = Sys.time()
        if(verbose) cat(paste('ts:',i,'- test series evaluation: done - ', endTime - startTime,'\n'))
                          
        resultTableAux
    }

    return(resultTable)
}
