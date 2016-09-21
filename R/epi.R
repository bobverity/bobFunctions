
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
#' SIS_deterministic
#'
#' Returns solution to deterministic SIS model using the \code{deSolve} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_deterministic <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {
	
	# load deSolve
	if (!"deSolve"%in%rownames(installed.packages()))
		install.packages('deSolve')
	require(deSolve)
	
	# set up parameters and initial conditions
	params <- c(beta=beta, r=r, N=N)
	state <- c(S=N-I_init, I=I_init)
	
	# define ode
	ode1 <- function(t, state, params) {
		with(as.list(c(state, params)), {
			# rate of change
			dS <- -beta*S*I/N + r*I
			dI <- beta*S*I/N - r*I
			
			# return the rate of change
			list(c(dS, dI))
		})
	}
	
	# solve ode
	output <- as.data.frame(ode(state, times, ode1, params))
	
	return(output)
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
#' SIR_deterministic
#'
#' Returns solution to deterministic SIR model using the \code{deSolve} package.
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
	
	# load deSolve
	if (!"deSolve"%in%rownames(installed.packages()))
		install.packages('deSolve')
	require(deSolve)
	
	# set up parameters and initial conditions
	params <- c(beta=beta, r=r, mu=mu, N=N)
	state <- c(S=N-I_init-R_init, I=I_init, R=R_init)
	
	# define ode
	ode1 <- function(t, state, params) {
		with(as.list(c(state, params)), {
			# rate of change
			dS <- -beta*S*I/N + mu*(I+R)
			dI <- beta*S*I/N - r*I - mu*I
			dR <- r*I - mu*R
			
			# return the rate of change
			list(c(dS, dI, dR))
		})
	}
	
	# solve ode
	output <- as.data.frame(ode(state, times, ode1, params))
	
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
#' SIR_stochastic_sync
#'
#' Draw from synchronous stochastic SIR model. Return state of the system at all time points at which an event occurs. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
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
#' SIR_delay_deterministic
#'
#' Returns solution to deterministic SIR model in which individuals are infectious for a fixed amount of time and there is no natural birth and death. Solves delay differential equation using the \code{deSolve} package. A set number of infectious individuals \code{I_init} are assumed to be seeded at time \code{t=0}, and will recover simultaneously at time \code{t=dur_inf}.
#'
#' @param beta contact rate.
#' @param dur_inf length of time in infectious state.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_delay_deterministic <- function(beta=0.5, dur_inf=4, I_init=10, R_init=0, N=1e3, times=0:100) {
	
	# load deSolve
	if (!"deSolve"%in%rownames(installed.packages()))
		install.packages('deSolve')
	require(deSolve)
	
	# set up parameters and initial conditions
	params <- c(beta=beta, dur_inf=dur_inf, N=N)
	state <- c(S=N-I_init-R_init, I=I_init, R=R_init)
	
	# define ode
	ode1 <- function(t, state, params) {
		with(as.list(c(state, params)), {
			# delay states
			if (t<dur_inf) {
				S_delay <- 0
				I_delay <- 0
			}
			else {
				S_delay <- lagvalue(t-dur_inf)[1]
				I_delay <- lagvalue(t-dur_inf)[2]
			}
			
			# rate of change
			dS <- -beta*S*I/N
			dI <- beta*S*I/N - beta*S_delay*I_delay/N
			dR <- beta*S_delay*I_delay/N
			
			# return the rate of change
			list(c(dS, dI, dR))
		})
	}
	
	# implement the recovery of initial infectious cohort with a discrete event
	event1 <- function(t, state, params) {
		with(as.list(state), {
			state[2] <- state[2] - I_init
			state[3] <- state[3] + I_init
			return(state)
		})
	}
	
	# solve ode
	output <- as.data.frame(dede(state, times, ode1, params, events=list(func=event1, time=dur_inf)))
	
	return(output)
}

# -----------------------------------
#' SLIR_deterministic
#'
#' Returns solution to deterministic SLIR model, where L is an incubation (lag) stage of defined length. Solves delay differential equation using the \code{deSolve} package.
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
	
	# load deSolve
	if (!"deSolve"%in%rownames(installed.packages()))
		install.packages('deSolve')
	require(deSolve)
	
	# set up parameters and initial conditions
	params <- c(beta=beta, dur_lag=dur_lag, N=N)
	state <- c(S=N-I_init-R_init, L=0, I=I_init, R=R_init)
	
	# define ode
	ode1 <- function(t, state, params) {
		with(as.list(c(state, params)), {
			# delay states
			if (t<dur_lag) {
				S_delay <- 0
				I_delay <- 0
			}
			else {
				S_delay <- lagvalue(t-dur_lag)[1]
				I_delay <- lagvalue(t-dur_lag)[3]
			}
			
			# rate of change
			dS <- -beta*S*I/N
			dL <- beta*S*I/N - beta*S_delay*I_delay/N
			dI <- beta*S_delay*I_delay/N - r*I
			dR <- r*I
			
			# return the rate of change
			list(c(dS, dL, dI, dR))
		})
	}
	
	# solve ode
	output <- as.data.frame(dede(state, times, ode1, params))
	
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