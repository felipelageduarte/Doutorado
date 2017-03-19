
#Include required Librarys
require(fNonlinear)
require(tseriesChaos)
require(nonlinearTseries)

seriesSize = 3000
dataFolder = 'data'

normalize <- function(values, a = -1, b = 1){
  max = max(values)
  min = min(values)
  return( a + (((values - min)*(b - a))/(max - min)) )
}

getDeterministicSeries <- function(comp, size = 1000){
  if(missing(comp)){
    stop("You must select one of those options for deterministic component: sine, lorenz, rossler, logistic or henon")
  }

  if(tolower(comp) == 'sine'){
    return(list( det.series = normalize(sin(2*pi*seq(0,9,len=size))),
                 det.sepDim = 5,
                 det.embDim = 2,
                 det.model  = "sine" ,
                 det.length = size)
           )
  } else if(tolower(comp) == 'lorenz') {
    return(
      list(det.series = normalize(lorentzSim(doplot = FALSE,
                            times  = seq(0, 50, by = (50/size)),
                            parms  = c(sigma = 16, r = 45.92, b = 4),
                            start  = c(-14, -13, 47)
      )[,2][1:size]),
      det.sepDim = 1,
      det.embDim = 3,
      det.model  = "lorenz",
      det.length = size)
    )
  } else if(tolower(comp) == 'rossler') {
      return(
        list(det.series = normalize(sim.cont(rossler.syst,
                            start.time=0,
                            end.time=650,
                            dt=650/size,
                            start.x=c(0,0,0),
                            parms=c(0.15, 0.2, 10))[1:size]),
             det.sepDim = 1,
             det.embDim = 3,
             det.model  = "rossler",
             det.length = size)
      )
  } else if(tolower(comp) == 'logistic') {
    return(
      list(det.series = normalize(logisticSim(n=size, parms = c(r = 3.8), start = 0.5, doplot = FALSE)),
           det.sepDim = 1,
           det.embDim = 2,
           det.model  = "logistic",
           det.length = size)
    )
  } else if(tolower(comp) == 'henon') {
    return(
      list(det.series = normalize(nonlinearTseries::henon(
                                start=c(-0.006423277,-0.473545134),
                                n.sample = size,
                                n.transient=10,
                                do.plot=FALSE)$x[1:size]),
           det.sepDim = 1,
           det.embDim = 2,
           det.model  = "henon",
           det.length = size)
    )
  } else {
    stop("You must select one of those options for deterministic component: sine, lorenz, rossler, logistic or henon")
  }
}

getStochasticSeries <- function(comp, params=list(), size = 1000){
  set.seed(7)
  if(missing(comp)){
    stop("You must select one of those options for stochastic component: zero, uniforme or normal")
  }

  if(tolower(comp) == 'zero'){
    return(list( sto.series = rep(0, size),
                 sto.model  = "zero" ,
                 sto.params = list(),
                 sto.length = size)
    )
  } else if(tolower(comp) == 'uniforme') {
    min = ifelse(is.null(params$min), -1, params$min)
    max = ifelse(is.null(params$max),  1, params$max)
    return(list( sto.series = runif(size, min, max),
                 sto.model  = paste("uniforme_",min,"_",max, sep=""),
                 sto.params = list(min=min, max=max),
                 sto.length = size)
    )
  } else if(tolower(comp) == 'normal') {
    mean = ifelse(is.null(params$mean), 0, params$mean)
    sd   = ifelse(is.null(params$sd), 1, params$sd)
    return(list( sto.series = rnorm(size, mean=mean, sd=sd),
                 sto.model  = paste("normal_",mean,"_",sd, sep=""),
                 sto.params = list(mean=mean, sd=sd),
                 sto.length = size)
    )
  } else {
    stop("You must select one of those options for stochastic component: zero, uniforme or normal")
  }
}      

timeSeriesFactor <- function(det.comp, sto.comp, sto.params=list(), size = 1000){
  det = getDeterministicSeries(det.comp, size)
  sto = getStochasticSeries(sto.comp, sto.params, size)
  tsObj = list(series = det$det.series + sto$sto.series, size = size)
  tsObj = c(tsObj, det, sto)
  return(tsObj)
}

set.seed(6426830)
det = list("sine", "lorenz", "rossler", "logistic", "henon")
sto = list(list(comp="zero"),
         list(comp="uniforme",params=list(min=-0.01, max=0.01)),
         list(comp="uniforme",params=list(min=-0.05, max=0.05)),
         list(comp="uniforme",params=list(min=-0.10, max=0.10)),
         list(comp="uniforme",params=list(min=-0.15, max=0.15)),
         list(comp="uniforme",params=list(min=-0.20, max=0.20)),
         list(comp="uniforme",params=list(min=-0.25, max=0.25)),
         list(comp="normal",params=list(sd=0.01)),
         list(comp="normal",params=list(sd=0.05)),
         list(comp="normal",params=list(sd=0.10)),
         list(comp="normal",params=list(sd=0.15)),
         list(comp="normal",params=list(sd=0.20)),
         list(comp="normal",params=list(sd=0.25))
         )

idx = 1
for(i in 1:length(det)){
    for(j in 1:length(sto)){
      seriesObj = timeSeriesFactor(det[[i]], sto[[j]]$comp, sto[[j]]$params, seriesSize)
      filename  = paste(dataFolder, sprintf('/series_%02d', idx),'.RData',sep='')
      save(seriesObj, file=filename)
      idx = idx + 1
    }
}

