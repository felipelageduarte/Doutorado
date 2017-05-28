
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

emdmiDec <- function(series, par){
  r.emd  = EMD::emd(series, 1:length(series), boundary = 'wave')
  phases = t(apply(r.emd$imf, 2, function(x){calcPhase(x)}))
  n.size = r.emd$nimf - 1

  # l = 1
  idx = c(1)
  if(n.size > 1){
    mi = c()
    for(i in 1:n.size){
      mi = c(mi, c3net::makemim(phases[i:(i+1),])[1,2])
    }
    idx = which(mi >= mean(mi))

    # mi = mapply(function(x,y) c3net::makemim(x,y),
    #             as.data.frame(phases[,1:n.size]),
    #             as.data.frame(phases[,2:(n.size+1)]))

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

