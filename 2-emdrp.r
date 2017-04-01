
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

waveletDec <- function(series, par){
  filter   = unlist(par[1])
  n.levels = as.numeric(unlist(par[2]))
  boundary = unlist(par[3])
  r.wavelet = wavelets::dwt(series,
                            filter = filter,
                            n.levels=n.levels,
                            boundary=boundary,
                            fast=TRUE)
  for (i in 1:length(r.wavelet@W)) {
    r.wavelet@W[[i]] = cbind(rep(0, length(r.wavelet@W[[i]])))
  }
  det = wavelets::idwt(r.wavelet)
  return(det)
}

params = expand.grid(
  detlevel = seq(0.01, 0.99, by=0.01),
  thresh   = seq(0.01, 0.99, by=0.01),
  delay    = 0:4,
  embedded = 0:4,
  stringsAsFactors = FALSE
)

resultTable = gridSearch(emdrpDec, params, seriesList, modelFolder, 'EMDRP', cores, TRUE)

write.csv(resultTable, file=paste(resultFolder,'/emdrp.csv', sep=''), row.names=FALSE)

