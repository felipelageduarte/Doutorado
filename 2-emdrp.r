
require(TSDecomposition)
require(foreach)
require(parallel)
require(doMC)
require(Rssa)
require(FNN)

source('utils.r')

cores = 1#detectCores(all.tests = FALSE, logical = TRUE)

dataFolder = 'data'
modelFolder  = 'model'
resultFolder = 'testResult'

seriesList = loadSeriesFile(dataFolder)

emdrpDec   <- function(series, params){
  detlevel = unlist(params[1])
  thresh   = unlist(params[2])
  delay    = unlist(params[3])
  embedded = unlist(params[4])
  emdrp    = rpemdDecomposition(series, detlevel, thresh, delay, embedded)
  return(emdrp@deterministic)
}

params = expand.grid(
  detlevel = seq(0.01, 0.99, by=0.05),
  thresh   = seq(0.01, 0.99, by=0.05),
  stringsAsFactors = FALSE
)

resultTable = gridSearch(emdrpDec, params, seriesList, modelFolder, 'emdrp', cores)

write.csv(resultTable, file=paste(resultFolder,'/emdrp.csv', sep=''), row.names=FALSE)

