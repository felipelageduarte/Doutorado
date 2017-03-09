md.dist <- function(d1, d2){
	if(ncol(d1) != ncol(d2)) stop("Arrays must have the same dimensions")
	if(nrow(d1) != nrow(d2)) stop("Arrays must have the same dimensions")

	d = d1 - d2
	sqrt(diag(d%*%t(d)))
}


mddl <-function(det.original, det.result){
	s1 = det.original
	s2 = det.result

	if(length(s1) > length(s2)){
		s.aux = s1
		s1 = s2
		s2 = s.aux
	}

	n.elem = length(s1)
	delta = length(s2) - length(s1)

	mddl = Inf
	for(i in 0:delta){
		d1     = s1
		d2     = s2[(1+i):(n.elem+i)]
		mddl = min(mddl, TSDecomposition::mddl(d1, d2, plot=FALSE))
	}
	return(mddl)
}


mda <- function(det.original, det.result, m, d){
	at1 = embedd(det.original, m=m, d=d)
	at2 = embedd(det.result, m=m, d=d)

	if(nrow(at1) > nrow(at2)){
		at.aux = at1
		at1 = at2
		at2 = at.aux
	}

	n.elem = nrow(at1)
	delta = nrow(at2) - nrow(at1)

	mda = Inf
	for(i in 0:delta){
		d1     = at1
		d2     = at2[(1+i):(n.elem+i),]
		p.dist = md.dist(d1,d2)
		mda    = min(mda, p.dist)
	}

	return(mda)
}

gridSearch <- function(F, params, seriesFile, tech.name, cores = 1){

	registerDoMC(cores)

	resultTable = foreach(i=1:length(seriesFile), .combine='rbind') %dopar% {
		s = seriesFile[i]
		resultTable.aux = matrix(ncol=6, nrow = 0)

		#load series object from RData file
		load(s)

		mddl.results = c()
		mda.results = c()

		for(p in 1:nrow(params)){

			result = c()

			#apply model on the time series
			tryCatch({
				result = F(seriesObj, params[p,])
			}, warning = function(w){
				print(w)
			}, error = function(e){
				print(e)
			}, finally = {
			})

			#---------------------------------------------------------------------------
			# This sequence of if-else is necessary because:
			# 1) Non all of parameter combination is valid and in this case the function 
			#    F will throw exception and return NULL as result. Instead of map all 
			#    valid combination we prefer to verify after code execution to keep it 
			#    more readable.
			# 2) We embedded the evaluation inside ForceDec execution in attempt to make 
			#    test more fast. In this situation the result will be a list with MDA 
			#    and MDDL values. 
			# 3) In all other non exception scenarios will return the deterministic 
			#    component and all evaluation must be done.
			#---------------------------------------------------------------------------
			if(is.null(result)){
				next
			} else if(is.list(result)){ 
				mddl.results = c(mddl.results, result$mddl)
				mda.results  = c(mda.results,  result$mda)

				for(k in 1:length(result$mddl)){
					par    = params[p,]
					par[2] = k

					resultTable.aux  = rbind(resultTable.aux, c( ((p-1)*length(result$mddl)) + k, 
																tech.name, 
																paste(par, collapse=","), 
																s, 
																result$mddl[k], 
																result$mda[k]))
				} 
			} else{
				mddl.r = mddl(seriesObj$det.series, result)
				mda.r  = mda (seriesObj$det.series, result, seriesObj$det.embDim, seriesObj$det.sepDim)
				mddl.results = c(mddl.results, mddl.r)
				mda.results  = c(mda.results,  mda.r)
				resultTable.aux  = rbind(resultTable.aux, c(p, tech.name, paste(params[p,],collapse=","), s, mddl.r, mda.r))
			}

			print(paste("Done:",i,'\n'))
		}

		#normalize results and calculate L2-norm
		l2n.results      = sqrt((mddl.results/max(mddl.results))^2 + (mda.results/max(mda.results))^2)
		best.result.idx  = which.min(l2n.results)
		best.result.idx  = ifelse(length(best.result.idx) == 0, 1, best.result.idx)

		model = list( model.name   = tech.name,
					 F            = F,
					 best.params  = as.list(strsplit(resultTable.aux[best.result.idx,3], ",")[[1]]),
					 mda          = mda.results[best.result.idx],
					 mddl         = mddl.results[best.result.idx]
					 )

		save(model, file=paste(modelFolder, '/', tech.name, '_', gsub("[^\\d]+", "", s, perl=TRUE) ,'.RData',sep=''))
		resultTable.aux
	}


	resultTable = data.frame(resultTable, stringsAsFactors = FALSE)
	colnames(resultTable) = c("id", "model", "model.params", "series", "mddl", "mda")
	resultTable = resultTable[order(resultTable$series, as.numeric(resultTable$id)),]
	resultTable$id = 1:nrow(resultTable)

	return(resultTable)
}

