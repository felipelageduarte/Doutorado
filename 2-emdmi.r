
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

emd.fixed <- function(xt, tt=NULL, tol=sd(xt)*0.1^2, max.sift=20, stoprule="type1",
                      boundary="periodic", sm="none", smlevels=c(1), spar=NULL, alpha=NULL,
                      check=FALSE, max.imf=10, plot.imf=FALSE, interm=NULL, weight=NULL){
  r.emd  = EMD::emd(xt, tt, tol, max.sift, stoprule, boundary, sm, smlevels,
                    spar, alpha, check, max.imf, plot.imf, interm, weight)
  r.emd$imf[which(abs(r.emd$imf) < 10^-12)] = 0
  r.emd$imf = r.emd$imf[,which(colSums(r.emd$imf) != 0)]
  r.emd$imf = data.frame(r.emd$imf)
  r.emd$nimf = ncol(r.emd$imf)
  r.emd$residue[which(r.emd$residue < (10^-12))] = 0
  r.emd
}

emdmiDec <- function(series, par){
  r.emd  = emd.fixed(xt = series,  tt = 1:length(series), boundary = 'wave')
  phases = apply(r.emd$imf, 2, function(x){calcPhase(x)})
  phases = cbind(phases, calcPhase(r.emd$residue))
  n.size = r.emd$nimf
  l = 1
  if(n.size > 1){
    mi = mapply(function(x,y) FNN::mutinfo(x,y),
                as.data.frame(phases[,1:n.size]),
                as.data.frame(phases[,2:(n.size+1)]))
    if(length(mi) > 1){
      miDiff = abs(diff(mi))
      l      = max(which(miDiff == max(miDiff)))+1
    }
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

