
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

ssaDec <- function(series, par){
  L        = as.numeric(unlist(par[1]))
  neig     = as.numeric(unlist(par[2]))
  kind     = unlist(par[3])

  #execute
  s = Rssa::ssa(series, L=L, neig=neig, kind=kind)
  r = reconstruct(s, groups = seq(1:(L/2)))

  #mutual information to separate deterministic components
  mi = c()
  for(i in 1:(length(r)-1))
    mi = c(mi, FNN::mutinfo(r[[i]], r[[i+1]]))

  det.idx = which.max(abs(diff(mi))) + 1

  #sum det. comp.
  detComp = rep(0, length(r[[1]]))
  for(i in 1:det.idx)
    detComp = detComp + r[[i]]

  return(detComp)
}

params = expand.grid(
  list(L = c(6:10,25,50),
  neig = c(3:10,25,50),
  kind = c('1d-ssa', 'toeplitz-ssa', 'mssa', '2d-ssa', 'nd-ssa')),
  stringsAsFactors = FALSE
)

resultTable = gridSearch(ssaDec, params, seriesList, modelFolder, 'SSA', cores)

write.csv(resultTable, file=paste(resultFolder,'/ssa.csv', sep=''), row.names=FALSE)



