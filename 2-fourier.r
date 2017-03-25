setwd("~/Dropbox/USP/Doutorado/ForceDecPaper")

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

fourierDec <- function(series, par){
    freq.cutoff = unlist(par[1])

    coeffs = fft(series)
    mags   = coeffs[1:(length(coeffs)/2)]
    mags   = 1+sqrt(Re(mags)^2+Im(mags)^2)
    if(freq.cutoff > length(mags))  stop('Frequence does not exist')
    o.idx  = order(mags, decreasing = T)
    idx    = (1:length(mags))[-o.idx[1:freq.cutoff]]

    coeffs[idx] = complex(real=0, imaginary=0)
    coeffs[length(coeffs) - idx + 1] = complex(real=0, imaginary=0)
    det = Re(fft(coeffs, inverse=T)) / length(series)

    return(det)
}

#all combination of possible parameters for Fourier Algorithm
params = expand.grid(
  cutoff = 1:450,
  stringsAsFactors = FALSE
)

resultTable = gridSearch(fourierDec, params, seriesList,
                         modelFolder, 'fourier', cores, TRUE)

write.csv(resultTable, file=paste(resultFolder,'/fourier.csv', sep=''), row.names=FALSE)



