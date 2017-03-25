
require(foreach)
require(parallel)
require(doMC)
require(FNN)

source('utils.r')

cores = detectCores(all.tests = FALSE, logical = TRUE)

seriesFolder = 'data'
modelFolder  = 'model'
resultFolder = 'testResult'

#select all series file in the data folder
seriesFile = list.files(path = seriesFolder, pattern = NULL, all.files = FALSE, full.names = TRUE, recursive = FALSE)

forceDec   <- function(seriesObj, par, test.execution=TRUE){
	k        = unlist(par[1])
	num.it   = unlist(par[2])
	epsilon  = unlist(par[3])
	delta    = unlist(par[4])
	delay    = unlist(seriesObj$det.sepDim + par[5])
	embedded = unlist(seriesObj$det.embDim + par[6])
	series   = unlist(seriesObj$series)

	mddl.results = c()
	mda.results = c()

	s.emb = embedd(seriesObj$series, m=embedded, d=delay)
	pos.i = s.emb
	nn = get.knn(s.emb, k=k+1)$nn.index #search for k-nearest neighbor
	nn = cbind(1:nrow(nn), nn) #place itself as a neighbor
	norm.simil = c()
	for (i in 1:nrow(nn)){
		neigh.pos  = s.emb[nn[i,],]
		neigh.dist = as.matrix(dist(neigh.pos))[1,]
		neigh.dist[which(neigh.dist < epsilon)] = 0
		similarity = 1 - neigh.dist/(max(neigh.dist)+epsilon)
		similarity[which(similarity < epsilon)] = 0
		norm.simil = rbind(norm.simil, similarity[-1]/(sum(similarity[-1])))
	}

	disp.vec = c()
	it = -1
	#this loop can be done parallel to accelerate processing
	for(it in 1:num.it){ 
		disp = 0
		for (i in 1:nrow(nn)){
			neigh.pos  = s.emb[nn[i,],]
			pos.i[i,]  = ((  delta   * neigh.pos[1,])
						  + ((1-delta)  * colSums(neigh.pos[-1,]*norm.simil[i,])))
			disp = disp + dist(rbind(pos.i[i,], s.emb[i,]))
		}
		disp = disp/nrow(nn)
		disp.vec = c(disp.vec, disp)
		#if(disp <= epsilon) break;
		s.emb = pos.i

		#--------------------------------------------------------------
		# This evaluation is running embedded in the forceDec function
		# in order to accelerate the test, indeed this code does not
		# belong to this function and must not be implemented in real
		# scenarios
		#--------------------------------------------------------------
		if(test.execution){
			mddl.r = mddl(seriesObj$det.series, s.emb[,1])
			mda.r  = mda (seriesObj$det.series, s.emb[,1], embedded, delay)

			mddl.results = c(mddl.results, mddl.r)
			mda.results  = c(mda.results,  mda.r)
		}
		#--------------------------------------------------------------
		# Evaluation code end.
		#--------------------------------------------------------------
	}

	return(list(mda = mda.results, mddl = mddl.results))
}

params = expand.grid(
					 k = c(2:10,15,20,25,30,25,40,45,50),
					 num.it = 30,
					 epsilon = 10^-6,
					 delta = c(0.01, seq(0.05, 1, by=0.05)),
					 delay = 0:4,
					 embedded = 0:4
					 )
nrow(params)

resultTable = gridSearch(forceDec, params, seriesFile, 'ForceDec', cores)

paste('FIM DOS TESTES!!!!!')

write.csv(resultTable, file=paste(resultFolder,'/forcedec.csv', sep=''))


