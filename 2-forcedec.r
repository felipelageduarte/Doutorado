
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

nn.similarity <- function(neigh.pos.i){
  neigh.dist     = as.matrix(dist(neigh.pos.i))[1,]
  similarity     = (1 - (neigh.dist/(max(neigh.dist)+10^-12)))
  similarity[-1] = similarity[-1]/sum(similarity[-1])
  return(similarity)
}

forceDec   <- function(series, par, test.execution=TRUE){
  k        = unlist(par[1])
  num.it   = unlist(par[2])
  epsilon  = unlist(par[3])
  delta    = unlist(par[4])
  delay    = unlist(par[5])
  embedded = unlist(par[6])

  s.emb = embedd(series, m=embedded, d=delay)
  nn = get.knn(s.emb, k=k)$nn.index #search for k-nearest neighbor
  nn = cbind(1:nrow(nn), nn) #place itself as a neighbor

  for(it in 1:num.it){
    neigh.pos  = lapply(split(nn, 1:nrow(nn)), function(x){s.emb[x,]})
    norm.simil = lapply(neigh.pos, function(x){nn.similarity(x)})
    pos = t(mapply(
      function(x,y, delta){
        ((delta * x[1,]) + ((1-delta)  * colSums(x[-1,]*y[-1])))
      },
      neigh.pos,
      norm.simil,
      MoreArgs = list(delta=delta)
    ))

    d = s.emb - pos
    disp = mean(sqrt(diag(d%*%t(d))))
    cat(paste(disp, '\n'))

    if(disp <= epsilon) break;
    s.emb = pos
  }
  return(s.emb[,1])
}

params = expand.grid(
  k = c(2:10,25,50),
  num.it = 30,
  epsilon = seq(0.001, 0.01,by=0.001),
  delta = c(0, 0.25, 0.5, 0.75, 1)
)

resultTable = gridSearch(waveletDec, params, seriesList, modelFolder, 'forcedec', cores)

write.csv(resultTable, file=paste(resultFolder,'/forcedec.csv', sep=''), row.names=FALSE)

