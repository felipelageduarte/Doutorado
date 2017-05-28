
source('utils.r')

cores = detectCores(all.tests = FALSE, logical = TRUE)

dataFolder = 'data'
modelFolder  = 'model'
resultFolder = 'testResult'

seriesList = loadSeriesFile(dataFolder)

dimSepMean <- function(s.emb, m, d){
  n = nrow(s.emb) + (m-1)*d
  idx = matrix(1:n, ncol=m, nrow=n)
  for(i in 2:m) idx[,i] = idx[,i] - (i-1)*d
  idx[which(idx < 0 | idx > nrow(s.emb))] = 0
  for(i in 1:nrow(idx)){
    s.emb[cbind(idx[i,], 1:m)] = mean(s.emb[cbind(idx[i,], 1:m)])
  }
  s.emb
}

nn.similarity <- function(neigh.pos.i){
  neigh.dist     = as.matrix(dist(neigh.pos.i))[1,]
  similarity     = (1 - (neigh.dist/(max(neigh.dist)+10^-12)))
  similarity[-1] = similarity[-1]/sum(similarity[-1])
  return(similarity)
}

forceDecCorrection   <- function(series, par){
  k        = as.numeric(unlist(par[1]))
  num.it   = as.numeric(unlist(par[2]))
  epsilon  = as.numeric(unlist(par[3]))
  delta    = as.numeric(unlist(par[4]))
  embedded = as.numeric(unlist(par[5]))
  delay    = as.numeric(unlist(par[6]))

  s.emb = tseriesChaos::embedd(series, m=embedded, d=delay)
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

    pos = dimSepMean(pos, embedded, delay)
    d = s.emb - pos
    disp = mean(sqrt(diag(d%*%t(d))))

    if(disp <= epsilon) break;
    s.emb = pos
  }

  return(toTimeSpace(s.emb, embedded, delay))
}

params = expand.grid(
  k = c(2:10,25,50),
  num.it = 30,
  epsilon = seq(0.001, 0.01,by=0.001),
  delta = c(0, 0.25, 0.5, 0.75, 1)
)

resultTable = gridSearch(forceDecCorrection, params, seriesList, modelFolder, 'forceDecCorrection', cores)

write.csv(resultTable, file=paste(resultFolder,'/forceDecCorrection.csv', sep=''), row.names=FALSE)

