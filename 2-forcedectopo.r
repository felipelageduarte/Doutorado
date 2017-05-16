
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

calcNewPosition <- function(x, y, delta){
  if(is.vector(x)) return(x)
  else if(is.vector(x[-1,])) return((delta * x[1,]) + ((1-delta) * x[-1,] * y[-1]))
  else return((delta * x[1,]) + ((1-delta) * colSums(x[-1,]*y[-1])))
}

forceDecTopo <- function(series, par){

  r        = as.numeric(unlist(par[1]))
  num.it   = as.numeric(unlist(par[2]))
  epsilon  = as.numeric(unlist(par[3]))
  delta    = as.numeric(unlist(par[4]))
  embedded = as.numeric(unlist(par[5]))
  delay    = as.numeric(unlist(par[6]))

  s.emb = tseriesChaos::embedd(series, m=embedded, d=delay)
  d.mat = as.matrix(dist(s.emb, diag = T, upper = T))
  #search for neighbor inside radio R
  nnr   = apply(d.mat, 1, function(x){as.integer(names(which(sort(x) < r)))})
  #hist(unlist(lapply(nnr, length)))

  for(it in 1:num.it){
    neigh.pos  = lapply(nnr, function(x){s.emb[x,]})
    norm.simil = lapply(neigh.pos, function(x){nn.similarity(x)})
    pos = t(mapply(calcNewPosition, neigh.pos, norm.simil,MoreArgs = list(delta=delta)))
    d = s.emb - pos
    disp = mean(sqrt(diag(d%*%t(d))))

    if(disp <= epsilon) break;
    s.emb = pos
  }

  return(toTimeSpace(s.emb, embedded, delay))
}

params = expand.grid(
  r = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5),
  num.it = 30,
  epsilon = seq(0.001, 0.01,by=0.001),
  delta = c(0, 0.25, 0.5, 0.75)
)

resultTable = gridSearch(forceDecTopo, params, seriesList, modelFolder, 'forcedectopo', cores)

write.csv(resultTable, file=paste(resultFolder,'/forcedectopo.csv', sep=''), row.names=FALSE)

