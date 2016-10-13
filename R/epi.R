
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

simQuantiles <- function(FUN="RM1_stochastic_sync", args=list(), reps=1e2, quantiles=c(0.05,0.5,0.95)) {
	
	# run function once to get dimensions and variable names
	testOutput <- do.call(FUN, args)
	varNames <- setdiff(names(testOutput),"time")
	
	# repeat simulation many times and store in array
	simArray <- array(0, dim=c(nrow(testOutput), length(varNames), reps))
	for (i in 1:reps) {
		simArray[,,i] <- as.matrix(do.call(FUN, args)[,varNames])
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
	output <- as.data.frame(suppressWarnings(dede(state, times, ode1, params, events=list(func=event1, time=dur_inf))))
	output <- output[match(times,output$time),]
	
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

# -----------------------------------
#' RM1_deterministic
#'
#' Returns solution to deterministic Ross-Macdonald model. Solves delay differential equation using the \code{deSolve} package.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param p mosquito probability of surviving one day.
#' @param g mosquito instantaneous death rate. g = -log(p) unless specified.
#' @param u intrinsic incubation period. The number of days from infection to infectiousness in a human host.
#' @param v extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param E_h initial number of infected but not infectious humans.
#' @param I_h initial number of infectious humans.
#' @param H human population size.
#' @param E_m initial number of infected but not yet infectious mosquitoes.
#' @param I_m initial number of infectious mosquitoes.
#' @param m ratio of adult female mosquitoes to humans. Population density of adult female mosquitoes is equal to M = m*H.
#' @param times vector of times at which output should be returned.
#'
#' @export

RM1_deterministic <- function(a=0.3, p=0.9, g=NULL, u=22, v=10, r=1/200, b=1, c=1, E_h=0, I_h=10, H=100, E_m=0, I_m=0, m=1, times=0:100) {
	
	# load deSolve
	require(deSolve)
	
	# calculate g from p if not specified
	if (is.null(g))
		g <- -log(p)
	
	# set up parameters and initial conditions
	params <- c(a=a, p=p, g=g, u=u, v=v, r=r, b=b, c=c, H=H, m=m)
	state <- c(S_h=H-E_h-I_h, E_h=E_h, I_h=I_h, S_m=m*H-E_m-I_m, E_m=E_m, I_m=I_m)
	
	# define ode
	ode1 <- function(t, state, params) {
		with(as.list(c(state, params)), {
			# delay states
			if (t<u) {	# u = human time from infection to infectious
				S_h_du <- 0
				I_m_du <- 0
			}
			else {
				S_h_du <- lagvalue(t-u)[1]
				I_m_du <- lagvalue(t-u)[6]
			}
			if (t<v) {	# v = mosquito time from infection to infectious
				I_h_dv <- 0
				S_m_dv <- 0
			}
			else {
				I_h_dv <- lagvalue(t-v)[3]
				S_m_dv <- lagvalue(t-v)[4]
			}
			
			# human rates of change
			dS_h <- -a*b*S_h*I_m/H + r*I_h
			dE_h <- a*b*S_h*I_m/H - a*b*S_h_du*I_m_du/H
			dI_h <- a*b*S_h_du*I_m_du/H - r*I_h
			
			# mosquito rates of change
			dS_m <- -a*c*S_m*I_h/H + g*(S_m+E_m+I_m) - g*S_m
			dE_m <- a*c*S_m*I_h/H - a*c*S_m_dv*I_h_dv/H*exp(-g*v) - g*E_m
			dI_m <- a*c*S_m_dv*I_h_dv/H*exp(-g*v) - g*I_m
			
			# return rates of change
			list(c(dS_h, dE_h, dI_h, dS_m, dE_m, dI_m))
		})
	}
	
	# implement initial progression from E_h and E_m states as discrete events
	df_events <- data.frame(
		var = c("E_h", "I_h", "E_m", "I_m"),
		time = c(u, u, v, v),
		value = c(-E_h, E_h, -E_m*exp(-g*v), E_m*exp(-g*v)),
		method = rep("add", 4)
		)
	
	# solve ode
	output <- as.data.frame(suppressWarnings(dede(state, times, ode1, params, events=list(data=df_events))))
	output <- output[match(times,output$time),]
	
	return(output)
}

# -----------------------------------
#' RM1_stochastic_async
#'
#' Draw from asynchronous stochastic Ross-Macdonald model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param p mosquito probability of surviving one day.
#' @param g mosquito instantaneous death rate. g = -log(p) unless specified.
#' @param u intrinsic incubation period. The number of days from infection to infectiousness in a human host.
#' @param v extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param E_h initial number of infected but not infectious humans.
#' @param I_h initial number of infectious humans.
#' @param H human population size.
#' @param E_m initial number of infected but not yet infectious mosquitoes.
#' @param I_m initial number of infectious mosquitoes.
#' @param m ratio of adult female mosquitoes to humans. Population density of adult female mosquitoes is equal to M = m*H.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

RM1_stochastic_async <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Em_init=0, Im_init=0, H=100, M=100, maxIterations=1e4) {
	
	# calculate g from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# run model
	args <- list(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Em_init, Im_init=Im_init, H=H, M=M, maxIterations=maxIterations)
	rawOutput <- RM1_stochastic_async_cpp(args)
	
	# format output object
	output <- as.data.frame(rawOutput)
	output$H <- H
	output$M <- M
	output <- subset(output, Sh>=0)
	return(output)
}

# -----------------------------------
#' RM1_stochastic_hybrid
#'
#' Draw from asynchronous stochastic Ross-Macdonald model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless specified.
#' @param u intrinsic incubation period. The number of days from infection to infectiousness in a human host.
#' @param v extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param Eh_init initial number of infected but not infectious humans.
#' @param Ih_init initial number of infectious humans.
#' @param Em_init initial number of infected but not yet infectious mosquitoes.
#' @param Im_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

RM1_stochastic_hybrid <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, H=100, Em_init=0, Im_init=0, M=100, times=0:100, maxIterations=1e4) {
	
	# calculate mu from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# run model
	args <- list(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Em_init, Im_init=Im_init, H=H, M=M, t_vec=times, maxIterations=maxIterations)
	rawOutput <- RM1_stochastic_hybrid_cpp(args)
	
	output <- cbind(time=times, as.data.frame(rawOutput))
	output$H <- H
	output$M <- M
	output[output<0] <- NA
	return(output)
}

# -----------------------------------
#' RM1_stochastic_sync
#'
#' Draw from synchronous stochastic Ross-Macdonald model. Return state of the system at known time points. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless specified.
#' @param u intrinsic incubation period. The number of days from infection to infectiousness in a human host.
#' @param v extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param Eh_init initial number of infected but not infectious humans.
#' @param Ih_init initial number of infectious humans.
#' @param Em_init initial number of infected but not yet infectious mosquitoes.
#' @param Im_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#'
#' @export

RM1_stochastic_sync <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Em_init=0, Im_init=0, H=100, M=100, times=0:100) {
	
	# calculate mu from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# check that all time intervals are the same
	delta_t <- unique(times[-1]-times[-length(times)])
	if ((max(delta_t)-min(delta_t))>1e-10)
		stop("all time intervals in the vector 'times' must be equal")
	delta_t <- times[2]-times[1]
	
	# check that lag times are greater than delta_t
	if (u<delta_t | v<delta_t)
		stop("time between steps must be greater than lag times (both u and v)")
	
	# run model
	args <- list(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Em_init, Im_init=Im_init, H=H, M=M, times=times)
	rawOutput <- RM1_stochastic_sync_cpp(args)
	
	output <- cbind(time=times, as.data.frame(rawOutput))
	output$H <- H
	output$M <- M
	output[output<0] <- NA
	return(output)
}

# -----------------------------------
#' RM2_deterministic
#'
#' Variation of RM1 model - in this version the mosquito population grows in logistic manner given growth and death rates and a carrying capacity. The carrying capacity can be varied at any given time point to allow for changing mosquito densities over time.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param mu mosquito instantaneous  death rate.
#' @param lambda mosquito instantaneous birth rate.
#' @param u intrinsic incubation period. The number of days from infection to infectiousness in a human host.
#' @param v extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param Eh_init initial number of infected but not infectious humans.
#' @param Ih_init initial number of infectious humans.
#' @param Em_init initial number of infected but not yet infectious mosquitoes.
#' @param Im_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M_init initial mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#' @param Ktimes vector of times at which carrying capacity is defined.
#' @param Kvalues vector of carrying capacities that come into action at \code{Ktimes}.
#'
#' @export

RM2_deterministic <- function(a=0.3, mu=0.1, lambda=0.2, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Em_init=0, Im_init=0, H=100, M_init=100, times=0:100, Ktimes=c(0,100), Kvalues=c(100,200)) {
	
	# load deSolve
	require(deSolve)
	
	# set up parameters and initial conditions
	params <- list(a=a, mu=mu, lambda=lambda, u=u, v=v, r=r, b=b, c=c, H=H, Ktimes=Ktimes, Kvalues=Kvalues, Sh_to_Eh=0, Eh_to_Ih=0, Ih_to_Sh=0, Sm_to_Em=0, Em_to_Im=0)
	state <- c(Sh=H-Eh_init-Ih_init, Eh=Eh_init, Ih=Ih_init, Sm=M_init-Em_init-Im_init, Em=Em_init, Im=Im_init, M=M_init)
	
	# get index positions of states
	Sh_index <- which(names(state)=="Sh")
	Ih_index <- which(names(state)=="Ih")
	Sm_index <- which(names(state)=="Sm")
	Im_index <- which(names(state)=="Im")
	
	# define ode
	ode1 <- function(t, state, params) {
		with(as.list(c(state, params)), {
			
			# delay states
			if (t<u) {	# u = human time from infection to infectious
				Sh_du <- 0
				Im_du <- 0
			}
			else {
				Sh_du <- lagvalue(t-u, Sh_index)
				Im_du <- lagvalue(t-u, Im_index)
			}
			if (t<v) {	# v = mosquito time from infection to infectious
				Ih_dv <- 0
				Sm_dv <- 0
			}
			else {
				Ih_dv <- lagvalue(t-v, Ih_index)
				Sm_dv <- lagvalue(t-v, Sm_index)
			}
			
			# find carrying capacity at this point in time
			w <- findInterval(t, Ktimes)
			w <- ifelse(w>0, w, 1)
			K <- Kvalues[w]
			
			# mosquito population rates of change
			dM <- lambda*M*(1-M/K) - mu*M
			
			# human flows between states
			Sh_to_Eh <- a*b*Sh*Im/H
			Eh_to_Ih <- a*b*Sh_du*Im_du/H
			Ih_to_Sh <- r*Ih
			
			# human rates of change			
			dSh <- -Sh_to_Eh + Ih_to_Sh
			dEh <- Sh_to_Eh - Eh_to_Ih
			dIh <- Eh_to_Ih - Ih_to_Sh
			
			# mosquito flows between states
			Sm_to_Em <- a*c*Sm*Ih/H
			Em_to_Im <- a*c*Sm_dv*Ih_dv/H*exp(-mu*v)
			
			# mosquito rates of change
			dSm <- -Sm_to_Em - mu*Sm + lambda*M*(1-M/K)
			dEm <- Sm_to_Em - Em_to_Im - mu*Em
			dIm <- Em_to_Im - mu*Im
			
			# return output
			list(c(dSh, dEh, dIh, dSm, dEm, dIm, dM), Sh_to_Eh=Sh_to_Eh, Eh_to_Ih=Eh_to_Ih, Ih_to_Sh=Ih_to_Sh, Sm_to_Em=Sm_to_Em, Em_to_Im=Em_to_Im, H=H, K=K)
		})
	}
	
	# implement initial progression from E_h and E_m states as discrete events
	df_events <- data.frame(
		var = c("Eh", "Ih", "Em", "Im"),
		time = c(u, u, v, v),
		value = c(-Eh_init, Eh_init, -Em_init*exp(-mu*v), Em_init*exp(-mu*v)),
		method = rep("add", 4)
		)
	
	# solve ode
	output <- as.data.frame(suppressWarnings(dede(state, times, ode1, params, events=list(data=df_events))))
	output <- subset(output, time%in%times)
	
	return(output)
}

# -----------------------------------
#' RM2_stochastic_sync
#'
#' Draw from a synchronous stochastic version of a particular Ross-Macdonald-style model (see \code{?RM2_deterministic} for details of the model). Return state of the system at known time points. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param mu mosquito instantaneous  death rate.
#' @param lambda mosquito instantaneous birth rate.
#' @param u intrinsic incubation period. The number of days from infection to infectiousness in a human host.
#' @param v extrinsic incubation period. The number of days from infection to infectiousness in a mosquito host.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param Eh_init initial number of infected but not infectious humans.
#' @param Ih_init initial number of infectious humans.
#' @param Em_init initial number of infected but not yet infectious mosquitoes.
#' @param Im_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M_init initial mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#' @param Ktimes vector of times at which carrying capacity is defined.
#' @param Kvalues vector of carrying capacities that come into action at \code{Ktimes}.
#'
#' @export
	
RM2_stochastic_sync <- function(a=0.3, mu=0.1, lambda=0.2, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Em_init=0, Im_init=0, H=100, M_init=100, times=0:100, Ktimes=c(0,100), Kvalues=c(100,200)) {
	
	# check that all time intervals are the same
	delta_t <- unique(times[-1]-times[-length(times)])
	if ((max(delta_t)-min(delta_t))>1e-10)
		stop("all time intervals in the vector 'times' must be equal")
	delta_t <- times[2]-times[1]
	
	# check that lag times are greater than delta_t
	if (u<delta_t | v<delta_t)
		stop("time between steps must be greater than lag times (both u and v)")
	
	# run model
	args <- list(a=a, mu=mu, lambda=lambda, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Em_init, Im_init=Im_init, H=H, M_init=M_init, times=times, Ktimes=Ktimes, Kvalues=Kvalues)
	rawOutput <-  RM2_stochastic_sync_cpp(args)
	
	output <- cbind(time=times, as.data.frame(rawOutput))
	output[output<0] <- NA
	return(output)
}
