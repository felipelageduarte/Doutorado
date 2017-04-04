
require(TSDecomposition)
require(foreach)
require(parallel)
require(doMC)
require(Rssa)
require(FNN)

source('utils.r')

cores = detectCores(all.tests = FALSE, logical = TRUE)

dataFolder = 'data'
modelFolder  = 'model'
resultFolder = 'testResult'

seriesList = loadSeriesFile(dataFolder)

calcPhase <- function(x){
  fft = fft(x)
  atan(Im(fft)/Re(fft))
}

emdmiDec2 <- function(series, par){
  r.emd  = emd(xt = series,  tt = 1:length(series), boundary = 'wave')
  phases = apply(r.emd$imf, 2, function(x){calcPhase(x)})
  phases = cbind(phases, calcPhase(r.emd$residue))
  n.size = r.emd$nimf
  mi     = mapply(function(x,y) FNN::mutinfo(x,y),
                  as.data.frame(phases[,1:n.size]),
                  as.data.frame(phases[,2:(n.size+1)]))
  l = 1
  if(length(mi) > 1){
    miDiff = abs(diff(mi))
    l      = max(which(miDiff == max(miDiff)))+1
  }
  idx    = l:r.emd$nimf
  if(length(idx) > 1){
    plot(rowSums(r.emd$imf[,idx]) + r.emd$residue)
    return(rowSums(r.emd$imf[,idx]) + r.emd$residue)
  } else{
    return(r.emd$imf[,idx] + r.emd$residue)
  }
}

resultTable = gridSearch(emdmiDec, NULL, seriesList, modelFolder, 'emdmi', cores, TRUE)

write.csv(resultTable, file=paste(resultFolder,'/wavelet.csv', sep=''), row.names=FALSE)

