md.dist <- function(d1, d2){
  d = d1 - d2
  return(sqrt(diag(d%*%t(d))))
}

mddl <-function(obs, pred){
  return(TSDecomposition::mddl(obs, pred, plot=FALSE))
}

mda <- function(obs, pred, m, d){
  at1 = tseriesChaos::embedd(pred, m=m, d=d)
  at2 = tseriesChaos::embedd(obs,  m=m, d=d)
  n   = nrow(at1)
  d   = nrow(at2) - nrow(at1)

  if(d == 0) return(mean(md.dist(at1, at2)))
  else return(min(lapply(0:d, function(x) mean(md.dist(at1, at2[(1+x):(n+x),])))))
}

mae <- function(obs, pred){
  n = length(pred)
  d = length(obs) - length(pred)
  if(d == 0)return(mean(abs(obs-pred)))
  else return(min(lapply(0:d,function(x) mean(abs(obs[(1+x):(n+x),]-pred)))))
}

rmse <- function(obs, pred){
  n = length(pred)
  d = length(obs) - length(pred)
  if(d == 0) return(sqrt(mean((obs-pred)^2)))
  else return(min(lapply(0:d,function(x) sqrt(mean((obs[(1+x):(n+x),]-pred)^2)))))
}

loadSeriesFile <- function(seriesFolder){
  #select all series file in the data folder
  seriesFile = list.files(path = seriesFolder, full.names = TRUE)
  seriesList = list()
  for(i in 1:length(seriesFile))
    seriesList[[i]] = get(load(seriesFile[i]))
  return(seriesList)
}

exec <- function(s, F, param){
  result = c()
  tryCatch({
      result = F(s, param)
  }, warning = function(w){
    write(paste('WARN:', w, '\n'), stderr())
  }, error = function(e){
    write(paste('ERRO:', e, '\n'), stderr())
  })
  if(is.null(result)) result = rep(0, length(s))
  return(result)
}

foreachParam <- function(s, F, params){
  return(t(apply(params, 1, function(x) exec(s, F, x))))
}

evaluate <- function(obs, pred, m, d){
  return(c(mddl(obs, pred),
           mda (obs, pred, m, d),
           mae (obs, pred),
           rmse(obs, pred)
  ))
}

evaluateResult <- function(obs, resultSeries, m, d){
	resultTable = matrix(ncol=5, nrow=nrow(resultSeries))
  resultTable[,1:4] = apply(resultSeries, 1, function(pred) evaluate(obs, pred, m, d))
  standMddl = (resultTable[,1] - mean(resultTable[,1]))/(sd(resultTable[,1]))
  standMda  = (resultTable[,2] - mean(resultTable[,2]))/(sd(resultTable[,2]))
  resultTable[,5] = sqrt(standMddl^2 + standMda^2)
  colnames(resultTable) = c('MDDL', 'MDA','MAE', 'RMSE', 'Dist')
  return(data.frame(resultTable))
}

gridSearch <- function(F, params, seriesList, modelFolder, techName,
                       cores = 1, v = FALSE){

  doMC::registerDoMC(cores)

  i=0
  resultTable =	foreach::foreach(i=1:length(seriesList), .combine='rbind') %dopar% {
    #load series object from RData file
    seriesObj = seriesList[[i]]
    m = unlist(seriesObj$det.embDim)
    d = unlist(seriesObj$det.sepDim)

    #Cross-validation out-of-time
    l = ceiling(0.7*seriesObj$size)
    trainIdx = 1:(l-1)
    testIdx  = l:seriesObj$size

    #Grid Search for better params
    if(v) cat(paste('ts:',i,'- params tunning: begin\n'))
    st = Sys.time()
    resultSeries = foreachParam(seriesObj$series[trainIdx], F, params)
    et = Sys.time()
    if(v) cat(paste('ts:',i,'- params tunning: done - ', et-st,'\n'))

    if(v) cat(paste('ts:',i,'- evaluation: begin','\n'))
    st = Sys.time()
    validTestIdx   = which(abs(rowSums(resultSeries)) > 0)
    resultTableAux = evaluateResult(seriesObj$det.series[trainIdx],
                                    resultSeries[validTestIdx,], m, d)
    bestParamIdx   = which.min(resultTableAux$Dist)
		resultTableAux = cbind(rep(techName,nrow(resultTableAux)),
                           apply(params[validTestIdx,], 1,
                                 function(x) paste(x, collapse=",")),
                           resultTableAux
                          )
		colnames(resultTableAux) = c('Tech','Param',  'MDDL', 'MDA','MAE', 'RMSE', 'Dist')
		resultTableAux = data.frame(resultTableAux)
    et = Sys.time()
    if(v) cat(paste('ts:',i,'- evaluation: done - ', et-st,'\n'))

    if(v) cat(paste('ts:',i,'- test series evaluation: begin','\n'))
    st = Sys.time()
    #run best model with test data and evaluate it
		bestP  = params[validTestIdx[bestParamIdx],]
		pred   = F(seriesObj$series[testIdx], bestP)
    obs    = seriesObj$det.series[testIdx]
    result = evaluate(obs, pred, seriesObj$det.embDim, seriesObj$det.sepDim)

    #save model result into model folder
    model = list( model.name = techName, F = F, best.param = bestP, eval = result )
    save(model, file=paste(modelFolder, '/', techName, '_', formatC(i, width=2, format='d', flag='0') ,'.RData',sep=''))
    et = Sys.time()
    if(v) cat(paste('ts:',i,'- test series evaluation: done - ', et-st,'\n'))
    resultTableAux
  }
  return(resultTable)
}
