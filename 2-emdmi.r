
source('utils.r')

cores = detectCores(all.tests = FALSE, logical = TRUE)

dataFolder = 'data'
modelFolder  = 'model'
resultFolder = 'testResult'

seriesList = loadSeriesFile(dataFolder)

calcPhase <- function(x){
  fft = fft(x)
  r   = atan(Im(fft)/Re(fft))
  r[which(is.nan(r))] = 0
  return(r)
}

emd.fixed <- function(series, tt){
  r.emd  = EMD::emd(series, tt, boundary = 'wave')
  r.emd$imf[which(abs(r.emd$imf) < 10^-12)] = 0
  r.emd$imf = r.emd$imf[,which(colSums(r.emd$imf) != 0)]
  r.emd$imf = data.frame(r.emd$imf)
  r.emd$nimf = ncol(r.emd$imf)
  r.emd$residue[which(r.emd$residue < (10^-12))] = 0
  r.emd
}

emdmiDec <- function(series, par){
  r.emd  = emd.fixed(series, 1:length(series))
  phases = apply(r.emd$imf, 2, function(x){calcPhase(x)})
  phases = cbind(phases, calcPhase(r.emd$residue))
  n.size = r.emd$nimf
  # l = 1
  idx = c(1)
  if(n.size > 1){
    mi = mapply(function(x,y) FNN::mutinfo(x,y),
                as.data.frame(phases[,1:n.size]),
                as.data.frame(phases[,2:(n.size+1)]))
    idx = which(mi >= mean(mi))
    # if(length(mi) > 1){
    #
    #   t = mean(mi)
    #   miDiff = abs(diff(mi))
    #   l      = max(which(miDiff == max(miDiff)))
    # }
  }
  # idx    = l:r.emd$nimf
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

