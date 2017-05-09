
source('utils.r')

cores = detectCores(all.tests = FALSE, logical = TRUE)

dataFolder = 'data'
modelFolder  = 'model'
resultFolder = 'testResult'

seriesList = loadSeriesFile(dataFolder)

calcNewPosition <- function(x1){
  cbind(x1[,1],apply(x1[,-1],2,function(x2){(smooth.spline(x2, df=(length(x2)/2)))$y}))
}

forceDecTimeDomain <- function(series, par){

  l        = as.numeric(unlist(par[1]))
  num.it   = as.numeric(unlist(par[2]))
  epsilon  = as.numeric(unlist(par[3]))
  embedded = as.numeric(unlist(par[4]))
  delay    = as.numeric(unlist(par[5]))

  s.emb = tseriesChaos::embedd(series, m=embedded, d=delay)
  idxs = lapply(1:nrow(s.emb), function(x) {
    idxs = ((x-l) : ((x+l)))
    return(idxs[which(idxs > 0 & idxs <= nrow(s.emb))])
  })

  for(it in 1:num.it){
    neigh.pos  = lapply(idxs, function(x){cbind(x,s.emb[x,])})
    spline.pos = lapply(neigh.pos, calcNewPosition)
    pos        = matrix(unlist(lapply(1:nrow(s.emb), function(x) {
      spline.pos[[x]][which(spline.pos[[x]][,1] == x),-1]
    })), ncol=embedded, byrow = T)
    d = s.emb - pos
    disp = mean(sqrt(diag(d%*%t(d))))

    if(disp <= epsilon) break;
    s.emb = pos
  }

  return(s.emb[,1])
}

params = expand.grid(
  l = c(4:10, 15, 25, 30),
  num.it = 30,
  epsilon = seq(0.001, 0.01,by=0.001)
)

resultTable = gridSearch(forceDecTimeDomain, params, seriesList, modelFolder, 'forcedectimedomain', cores)

write.csv(resultTable, file=paste(resultFolder,'/forcedectimedomain.csv', sep=''), row.names=FALSE)

