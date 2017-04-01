
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

#all combination of possible parameters for wavelets Algorithm
hyperparameters = expand.grid(
  filters = c("haar", "d4", "d6", "d8", "d10", "d12", "d14", "d16", "d18", "d20",#Daubechies
              "la8", "la10", "la12", "la14", "la16", "la18", "la20", #Least Asymetric
              "bl14", "bl18", "bl20", #Best Localized
              "c6", "c12", "c18", "c24", "c30"), #Coiflet
  n.levels = 1:50,
  boundarys = c("periodic","reflection"),
  stringsAsFactors = FALSE
)

resultTable = gridSearch(waveletDec, hyperparameters, seriesList, modelFolder, 'wavelet', cores, TRUE)

write.csv(resultTable, file=paste(resultFolder,'/wavelet.csv', sep=''), row.names=FALSE)

