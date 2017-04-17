
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

emdmiDec <- function(series, par){
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
    return(rowSums(r.emd$imf[,idx]) + r.emd$residue)
  } else{
    return(r.emd$imf[,idx] + r.emd$residue)
  }
}

params = expand.grid(
  none = '',
  stringsAsFactors = FALSE
)

resultTable = gridSearch(emdmiDec, params, seriesList, modelFolder, 'emdmi', cores)

write.csv(resultTable, file=paste(resultFolder,'/emdmi.csv', sep=''), row.names=FALSE)

