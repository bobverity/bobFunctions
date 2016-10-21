
# -----------------------------------
#' safeRead
#'
#' Reads in text file (tab-delimited by default) and replaces bad characters with replacements.
#' \cr\cr
#' WARNING - can overwrite original file with new values if \code{overwrite==TRUE}, although this is not the default.
#'
#' @export

safeRead <- function(fileName, delim="\t", useHeader=TRUE, report=TRUE, badCharacters=c("-1.#IND","1.#INF","-999"), replacements=c(0,0,0), overwrite=FALSE) {
	
	# read in raw data. All values are read in as characters at this stage, and headers are ignored
	data <- read.delim(fileName, sep=delim, header=F, colClasses="character")
	
	# replace bad characters
	badCount <- 0
	for (j in 1:length(badCharacters)) {
		badCount <- badCount + sum(data==badCharacters[j])
		data[data==badCharacters[j]] <- replacements[j]
	}
	
	# report to console
	if (report & badCount>0) {
		cat(paste(fileName,": ",badCount," bad characters replaced\n",sep=""))
	}
	
	# now convert variables to sensible types 
	if (useHeader) {
		df <- data.frame(2:nrow(data))
		for (i in 1:ncol(data)) {
			df <- cbind(df, type.convert(data[-1,i]))
		}
		df <- df[,-1]
		names(df) <- as.character(data[1,])
	} else {
		df <- data.frame(1:nrow(data))
		for (i in 1:ncol(data)) {
			df <- cbind(df, type.convert(data[,i]))
		}
		df <- df[,-1]
		names(df) <- paste("X",1:ncol(data),sep="")
	}
	
	# return or write to file
	if (overwrite) {
		write.table(df, fileName, sep=delim, row.names=FALSE, col.names=useHeader, quote=FALSE)
	} else {
		return(df)
	}
}

# -----------------------------------
#' merge.SpatialPolygonsDataFrame
#'
#' Reads in a shape file and data frame to be merged with the data in this shapefile. Merges while preserving the order of objects (ordinary merge operation causes polygons to become disassociated with data).
#'
#' @export

merge.SpatialPolygonsDataFrame <- function(shp, df) {

	# load sp (I think this is the correct package?)
	if (!"sp"%in%rownames(installed.packages()))
		install.packages('sp')
	require(sp)

	# extract data from shapefile and add key
	shp_data <- shp@data
	shp_data$mergeKey <- 1:nrow(shp_data)
	
	# merge with df and sort based on key
	m <- merge(shp_data,df,all.x=T)
	m <- m[order(m$mergeKey),]
	m <- subset(m,select=-mergeKey)
	
	# fix row names
	row.names(m) <- row.names(shp)
	
	# make final SpatialPolygonsDataFrame object
	s <- SpatialPolygonsDataFrame(geometry(shp), m)
	return(s)
}

# -----------------------------------
#' getPolyArea
#'
#' Reads in a shape file and extracts area of every polygon.
#'
#' @export

getPolyArea <- function(shp) {

	polys <- slot(shp, "polygons")
	output <- rep(NA,length(polys))
	for (i in 1:length(polys)) {
		output[i] <- slot(polys[[i]], "area")
	}
	return(output)
}

# -----------------------------------
#' rateRatio
#'
#' Computes point estimate and upper and lower confidence intervals on a ratio of rates. Default method is to enter raw counts and time periods, but if entering rates simply set time1=1 and time2=1.
#'
#' @export

rateRatio <- function(count1, time1, count2, time2, alpha=0.05) {
	
	point <- time2/time1*count1/count2
	if (count1==0 | time1==0 | count2==0 | time2==0) {
		if (count1==0) {
			return(list(point=0,LL=NaN,UL=NaN))
		} else {
			return(list(point=NaN,LL=NaN,UL=NaN))
		}
	}
	LL <- time2/time1*count1/(count2+1)*qf(0.025,2*(count2+1),2*count1)
	UL <- time2/time1*(count1+1)/count2*qf(1-0.025,2*(count1+1),2*count2)
	
	return(list(point=point,LL=LL,UL=UL))
}

# -----------------------------------
#' simQuantiles
#'
#' Runs a given stochastic simulation function many times, computing the mean and quantiles over replicates. Note that this method will only work with simulations that have a fixed time step, i.e. synchronous or hybrid simulations, and not with asynchronous simulations. In the hybrid case the maxIterations limit cannot be reached in any simulation.
#'
#' @param FUN the stochastic simulation function to use.
#' @param args a list of arguments to the function.
#' @param reps number of times to repeat the stochastic simulation.
#' @param quantiles which quantiles to compute over replicates.
#'
#' @export

simQuantiles <- function(FUN="SIS_stochastic_hybrid", args=list(), reps=1e2, quantiles=c(0.05,0.5,0.95)) {
	
	# run function once to get dimensions and variable names
	testOutput <- do.call(FUN, args)
	varNames <- setdiff(names(testOutput),"time")
	
	# repeat simulation many times and store in array
	simArray <- array(0, dim=c(nrow(testOutput), length(varNames), reps))
	simArray[,,1] <- as.matrix(testOutput[,varNames])
	if (reps>1) {
		for (i in 2:reps) {
			simArray[,,i] <- as.matrix(do.call(FUN, args)[,varNames])
		}
	}
	
	# compute mean and quantiles over replicates and store in data frame
	df <- data.frame(time=testOutput$time)
	for (i in 1:length(varNames)) {
		m <- rowMeans(simArray[,i,,drop=FALSE])
		q <- apply(simArray[,i,,drop=FALSE], 1, quantile, prob=quantiles)
		df_new <- as.data.frame(t(rbind(m,q)))
		names(df_new) <- paste(varNames[i], c("mean", paste("Q", quantiles, sep="")), sep="_")
		df <- cbind(df, df_new)
	}
	
	# return summary data frame
	return(df)
}

# -----------------------------------
#' SIS_analytical
#'
#' Returns analytical solution to deterministic SIS model. At equilibrium the number of infectives is given by \eqn{I* = N(1 - r/beta)}.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_analytical <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {
	
	I <- (beta-r)/(beta/N + (beta-r-beta*I_init/N)/I_init*exp(-(beta-r)*times))
	
	output <- data.frame(time=times, S=N-I, I=I)
	return(output)
}

# -----------------------------------
#' SIS_deterministic
#'
#' Returns solution to deterministic SIS model using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_deterministic <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {

	# solve ode	
	mod <- SIS_deterministic_odin(beta=beta, r=r, I_init=I_init, N=N)
	output <- as.data.frame(mod$run(times))
	names(output)[1] <- 'time'
	
	return(output)
}

# -----------------------------------
#' SIS_stochastic_async
#'
#' Draw from asynchronous stochastic SIS model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIS_stochastic_async <- function(beta=1, r=0.25, I_init=100, N=1e3, maxIterations=1e4) {
	
	# run model
	args <- list(beta=beta, r=r, I_start=I_init, N=N, maxIterations=maxIterations)
	rawOutput <- SIS_stochastic_async_cpp(args)
	
	# format output object
	I <- rawOutput$I
	S <- N-I
	t_vec <- rawOutput$t
	output <- data.frame(time=t_vec,S=S,I=I)
	output <- subset(output, I>=0)

	return(output)
}

# -----------------------------------
#' SIS_stochastic_hybrid
#'
#' Draw from stochastic SIS model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIS_stochastic_hybrid <- function(beta=1, r=0.25, I_init=100, N=1e3, times=0:100, maxIterations=1e4) {
	
	# run model
	args <- list(beta=beta, r=r, I_start=I_init, N=N, t_vec=times, maxIterations=maxIterations)
	rawOutput <- SIS_stochastic_hybrid_cpp(args)
	
	# format output object
	I <- rawOutput$I
	S <- N-I
	S[I<0] <- NA
	I[I<0] <- NA
	output <- data.frame(time=times,S=S,I=I)

	return(output)
}

# -----------------------------------
#' SIS_stochastic_sync
#'
#' Draw from synchronous stochastic SIS model. Return infectives at known time points. Note that results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_stochastic_sync <- function(beta=1, r=0.25, I_init=100, N=1e3, times=0:100) {
	
	# run model
	args <- list(beta=beta, r=r, I_start=I_init, N=N, t_vec=times)
	rawOutput <- SIS_stochastic_sync_cpp(args)
	
	# format output object
	I <- rawOutput$I
	S <- N-I
	S[I<0] <- NA
	I[I<0] <- NA
	output <- data.frame(time=times,S=S,I=I)

	return(output)
}

# -----------------------------------
#' SIR_deterministic
#'
#' Returns solution to deterministic SIR model using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_deterministic <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100) {

	# solve ode	
	mod <- SIR_deterministic_odin(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N)
	output <- as.data.frame(mod$run(times))
	names(output)[1] <- 'time'
	
	return(output)
}

# -----------------------------------
#' SIR_stochastic_async
#'
#' Draw from asynchronous stochastic SIR model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIR_stochastic_async <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, maxIterations=1e4) {
	
	# run model
	args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, maxIterations=maxIterations)
	rawOutput <- SIR_stochastic_async_cpp(args)
	
	# format output object
	t_vec <- rawOutput$t
	S <- rawOutput$S
	I <- rawOutput$I
	R <- rawOutput$R
	output <- data.frame(time=t_vec, S=S, I=I, R=R)
	output <- subset(output, I>=0)

	return(output)
}

# -----------------------------------
#' SIR_stochastic_hybrid
#'
#' Draw from stochastic SIR model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIR_stochastic_hybrid <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100, maxIterations=1e4) {
	
	# run model
	args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, t_vec=times, maxIterations=maxIterations)
	rawOutput <- SIR_stochastic_hybrid_cpp(args)
	
	# format output object
	S <- rawOutput$S
	I <- rawOutput$I
	R <- rawOutput$R
	S[S<0] <- NA
	I[I<0] <- NA
	R[R<0] <- NA
	output <- data.frame(time=times, S=S, I=I, R=R)

	return(output)
}

# -----------------------------------
#' SIR_stochastic_sync
#'
#' Draw from synchronous stochastic SIR model. Return state of the system at known time points. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_stochastic_sync <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100) {
	
	# run model
	args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, t_vec=times)
	rawOutput <- SIR_stochastic_sync_cpp(args)
	
	# format output object
	S <- rawOutput$S
	I <- rawOutput$I
	R <- rawOutput$R
	S[S<0] <- NA
	I[I<0] <- NA
	R[R<0] <- NA
	output <- data.frame(time=times, S=S, I=I, R=R)

	return(output)
}

# -----------------------------------
#' SLIR_deterministic
#'
#' Returns solution to deterministic SLIR model, where L is an incubation (lag) stage of defined length. Solves delay differential equation using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SLIR_deterministic <- function(beta=0.5, dur_lag=1, r=0.25, I_init=10, R_init=0, N=1e3, times=0:100) {
	
	# solve ode
	mod <- SLIR_deterministic_odin(beta=beta, dur_lag=dur_lag, r=0.25, I_init=I_init, R_init=R_init, N=N)
	output <- as.data.frame(mod$run(times))
	names(output)[1] <- 'time'
	
	return(output)
}

# -----------------------------------
#' SLIR_stochastic_async
#'
#' Draw from asynchronous stochastic SLIR model, where L is an incubation (lag) stage of defined length. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SLIR_stochastic_async <- function(beta=1, dur_lag=1, r=0.25, I_init=100, R_init=0, N=1e3, maxIterations=1e4) {
	
	# run model
	args <- list(beta=beta, dur_lag=dur_lag, r=r, I_init=I_init, R_init=R_init, N=N, maxIterations=maxIterations)
	rawOutput <- SLIR_stochastic_async_cpp(args)
	
	# format output object
	t_vec <- rawOutput$t
	S <- rawOutput$S
	L <- rawOutput$L
	I <- rawOutput$I
	R <- rawOutput$R
	output <- data.frame(time=t_vec, S=S, L=L, I=I, R=R)
	output <- subset(output, I>=0)

	return(output)
}

# -----------------------------------
#' SLIR_stochastic_hybrid
#'
#' Draw from stochastic SLIR model, where L is an incubation (lag) stage of defined length, using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SLIR_stochastic_hybrid <- function(beta=1, dur_lag=1, r=0.25, I_init=100, R_init=0, N=1e3, times=0:100, maxIterations=1e4) {
	
	# run model
	args <- list(beta=beta, dur_lag=dur_lag, r=r, I_init=I_init, R_init=R_init, N=N, t_vec=times, maxIterations=maxIterations)
	rawOutput <- SLIR_stochastic_hybrid_cpp(args)
	
	# format output object
	S <- rawOutput$S
	L <- rawOutput$L
	I <- rawOutput$I
	R <- rawOutput$R
	S[S<0] <- NA
	L[L<0] <- NA
	I[I<0] <- NA
	R[R<0] <- NA
	output <- data.frame(time=times, S=S, L=L, I=I, R=R)

	return(output)
}
