
require(nonlinearTseries)
require(TSDecomposition)
require(tseriesChaos)
require(fNonlinear)
require(parallel)
require(foreach)
require(doMC)
require(Rssa)
require(FNN)

standardize <- function(x){
  if(length(x) == 1) return(x)
  else return((x-mean(x))/(sd(x)))
}
normalize <- function(values, a = -1, b = 1){
  max = max(values)
  min = min(values)
  return( a + (((values - min)*(b - a))/(max - min)) )
}

md.dist <- function(d1, d2){
  d = d1 - d2
  return(sqrt(mean(diag(d%*%t(d)))))
}
mddl <-function(obs, pred){
  return(rnorm(1))
  dtw = dtw(obs, pred)
  return(distanceToDiagonal(dtw$index1, dtw$index2, length(obs)))
}
mae.md <- function(obs, pred, m, d){
  at1 = tseriesChaos::embedd(pred, m=m, d=d)
  at2 = tseriesChaos::embedd(obs,  m=m, d=d)
  n   = nrow(at1)
  d   = nrow(at2) - nrow(at1)

  if(d == 0) return(mean(abs(at1 - at2)))
  else return(min(lapply(0:d, function(x) mean(abs(at1 - at2[(1+x):(n+x),])))))
}
rmse.md <- function(obs, pred, m, d){
  at1 = tseriesChaos::embedd(pred, m=m, d=d)
  at2 = tseriesChaos::embedd(obs,  m=m, d=d)
  n   = nrow(at1)
  d   = nrow(at2) - nrow(at1)

  if(d == 0) return(mean(md.dist(at1, at2)))
  else return(min(lapply(0:d, function(x) md.dist(at1, at2[(1+x):(n+x),]))))
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
evaluateMetrics <- function(obs, pred, m, d){
  return(c(mddl(obs, pred),
           mae.md(obs, pred, m, d),
           rmse.md(obs, pred, m, d),
           mae(obs, pred),
           rmse(obs, pred)
  ))
}
evaluateResult <- function(obs, resultSeries, params, techName, testId){
  resultTable = data.frame(testId   = numeric(0),  tech  = character(0),
                           paramIdx = numeric(0), param  = character(0),
                           mddl     = numeric(0), mae_md = numeric(0),
                           rmse_md  = numeric(0), mae    = numeric(0),
                           rmse     = numeric(0), dist   = numeric(0))

  validTestIdx = which(abs(rowSums(resultSeries)) > 0)
  resultSeries = matrix(resultSeries[validTestIdx,], ncol=ncol(resultSeries))

  m = unique(params$m)
  d = unique(params$d)

  er = apply(resultSeries, 1, function(pred) evaluateMetrics(obs, pred, m, d))
  resultTable[1:length(validTestIdx),5:9] = er
  resultTable$dist = sqrt(standardize(resultTable$mddl)^2 +
                          standardize(resultTable$rmse_md)^2)
  resultTable$testId = testId
  resultTable$tech = techName
  resultTable$paramIdx = validTestIdx
  resultTable$param = (apply(format(params), 1, paste, collapse=","))[validTestIdx]

  return(resultTable)
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
  if(is.null(params)) return(exec(s, F, NULL))
  return(t(apply(params, 1, function(x) exec(s, F, x))))
}

gridSearch <- function(F, params, seriesList, modelFolder, techName, cores = 1){
  resultTable =	foreach::foreach(i=1:3#length(seriesList)
                                 , .combine='rbind') %dopar% {
    st = Sys.time()
    cat(paste('ts:',i,'- begin\n'))

    #load series object from RData file
    seriesObj = seriesList[[i]]
    params$m = unlist(seriesObj$det.embDim)
    params$d = unlist(seriesObj$det.sepDim)

    #Grid Search for better params
    resultSeries = foreachParam(seriesObj$series, F, params)

    #evaluate results with know deterministic component
    det.comp = seriesObj$det.series
    rTable   = evaluateResult(det.comp, resultSeries, params, techName, i)
    bestIdx  = which.min(rTable$dist)

    #save model result into model folder
    model = list( model.name = techName, F = F,
                  best.param = params[rTable$ParamIdx[bestIdx],],
                  eval = rTable[bestIdx,] )
    save(model, file=paste(modelFolder, '/', techName, '_',
                           formatC(i, width=2, format='d', flag='0') ,'.RData',sep=''))

    et = Sys.time()
    cat(paste('ts:',i,'- done - ', et-st,'\n'))
    rTable
  }
  return(resultTable)
}
