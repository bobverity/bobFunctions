
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
SIS_deterministic_odin <- odin::odin({
	
	# derivatives
	deriv(S) <- -beta*S*I/N + r*I
	deriv(I) <- beta*S*I/N - r*I
	
	# initial conditions
	initial(S) <- N - I_init
	initial(I) <- I_init
	
	# parameters
	beta <- user()
	r <- user()
	I_init <- user()
	N <- user()
})

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
SIR_deterministic_odin <- odin::odin({
	
	# derivatives
	deriv(S) <- -beta*S*I/N + mu*(S+I+R) - mu*S
	deriv(I) <- beta*S*I/N - r*I - mu*I
	deriv(R) <- r*I - mu*R
	
	# initial conditions
	initial(S) <- N - I_init - R_init
	initial(I) <- I_init
	initial(R) <- R_init
	
	# parameters
	beta <- user()
	r <- user()
	mu <- user()
	I_init <- user()
	R_init <- user()
	N <- user()
})

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
#' Returns solution to deterministic SIR model in which individuals are infectious for a fixed amount of time and there is no natural birth and death. Solves delay differential equation using the \code{odin} package. A set number of infectious individuals \code{I_init} are assumed to be seeded at time \code{t=0}, and will recover simultaneously at time \code{t=dur_inf}.
#'
#' @param beta contact rate.
#' @param dur_inf length of time in infectious state.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_delay_deterministic2 <- function(beta=0.5, dur_inf=4, I_init=10, R_init=0, N=1e3, times=0:100) {

	# solve ode
	mod <- SIR_delay_deterministic_odin(beta=beta, dur_inf=dur_inf, I_init=I_init, R_init=R_init, N=N)
	output <- as.data.frame(mod$run(times))
	names(output)[1] <- 'time'
	
	return(output)
}
SIR_delay_deterministic_odin <- odin::odin({
	
	# delay states
	#S_delay <- delay(S, dur_inf, N - I_init - R_init)
	#I_delay <- delay(I, dur_inf, I_init)
	
	S_delay <- delay(S, dur_inf, 0)
	I_delay <- delay(I, dur_inf, 0)
	
	# derivatives
	deriv(S) <- -beta*S*I/N
	deriv(I) <- beta*S*I/N - beta*S_delay*I_delay/N
	deriv(R) <- beta*S_delay*I_delay/N
	
	# initial conditions
	initial(S) <- N - I_init - R_init
	initial(I) <- I_init
	initial(R) <- R_init
	
	# parameters
	beta <- user()
	dur_inf <- user()
	I_init <- user()
	R_init <- user()
	N <- user()
})

SIR_delay_deterministic <- function(beta=0.5, dur_inf=4, I_init=10, R_init=0, N=1e3, times=0:100) {
	
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
SLIR_deterministic_odin <- odin::odin({
	
	# delay states
	S_delay <- delay(S, dur_lag, 0)
	I_delay <- delay(I, dur_lag, 0)
	
	# derivatives
	deriv(S) <- -beta*S*I/N
	deriv(L) <- beta*S*I/N - beta*S_delay*I_delay/N
	deriv(I) <- beta*S_delay*I_delay/N - r*I
	deriv(R) <- r*I
	
	# initial conditions
	initial(S) <- N - I_init - R_init
	initial(L) <- 0
	initial(I) <- I_init
	initial(R) <- R_init
	
	# parameters
	beta <- user()
	dur_lag <- user()
	r <- user()
	I_init <- user()
	R_init <- user()
	N <- user()
})

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
#' Returns solution to deterministic Ross-Macdonald model. Solves delay differential equation using the \code{odin} package.
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
#' @param Ev_init initial number of infected but not yet infectious mosquitoes.
#' @param Iv_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#'
#' @export

RM1_deterministic <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Ev_init=0, Iv_init=0, H=100, M=100, times=0:100) {
	
	# calculate mu from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# solve ode
	mod <- RM1_deterministic_odin(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Ev_init=Ev_init, Iv_init=Iv_init, H=H, M=M)
	output <- as.data.frame(mod$run(times))
	names(output)[1] <- 'time'
	
	return(output)
}
RM1_deterministic_odin <- odin::odin({
	
	# delay states
	Sh_du <- delay(Sh, u, 0)
	Iv_du <- delay(Iv, u, 0)
	Sv_dv <- delay(Sv, v, 0)
	Ih_dv <- delay(Ih, v, 0)
	
	# human flows between states
	EIR <- a*Iv/H
	FOI <- b*EIR
	
	Sh_to_Eh <- FOI*Sh
	Eh_to_Ih <- (a*b*Iv_du/H)*Sh_du # equivalent to delayed force of infection * delayed Sh
	Ih_to_Sh <- r*Ih
	
	# human rates of change
	deriv(Sh) <- -Sh_to_Eh + Ih_to_Sh
	deriv(Eh) <- Sh_to_Eh - Eh_to_Ih
	deriv(Ih) <- Eh_to_Ih - Ih_to_Sh
	
	# mosquito flows between states
	Sv_to_Ev <- a*c*Sv*Ih/H
	Ev_to_Iv <- a*c*Sv_dv*Ih_dv/H*exp(-mu*v)
	
	# mosquito rates of change
	deriv(Sv) <- -Sv_to_Ev + mu*M - mu*Sv
	deriv(Ev) <- Sv_to_Ev - Ev_to_Iv - mu*Ev
	deriv(Iv) <- Ev_to_Iv - mu*Iv
	
	# initial conditions
	initial(Sh) <- H - Eh_init - Ih_init
	initial(Eh) <- Eh_init
	initial(Ih) <- Ih_init
	
	initial(Sv) <- M - Ev_init - Iv_init
	initial(Ev) <- Ev_init
	initial(Iv) <- Iv_init
	
	# parameters
	a <- user()
	mu <- user()
	u <- user()
	v <- user()
	r <- user()
	b <- user()
	c <- user()
	Eh_init <- user()
	Ih_init <- user()
	Ev_init <- user()
	Iv_init <- user()
	H <- user()
	M <- user()
	
	# other outputs
	output(H) <- H
	output(M) <- M
	output(EIR) <- EIR
	output(Sh_to_Eh) <- Sh_to_Eh
	output(Eh_to_Ih) <- Eh_to_Ih
	output(Ih_to_Sh) <- Ih_to_Sh
	output(Sv_to_Ev) <- Sv_to_Ev
	output(Ev_to_Iv) <- Ev_to_Iv
})

# -----------------------------------
#' RM1_stochastic_async
#'
#' Draw from asynchronous stochastic Ross-Macdonald model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
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
#' @param Ev_init initial number of infected but not yet infectious mosquitoes.
#' @param Iv_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

RM1_stochastic_async <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Ev_init=0, Iv_init=0, H=100, M=100, maxIterations=1e4) {
	
	# calculate g from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# run model
	args <- list(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Ev_init, Im_init=Iv_init, H=H, M=M, maxIterations=maxIterations)
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
#' @param Ev_init initial number of infected but not yet infectious mosquitoes.
#' @param Iv_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

RM1_stochastic_hybrid <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Ev_init=0, Iv_init=0, H=100, M=100, times=0:100, maxIterations=1e4) {
	
	# calculate mu from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# run model
	args <- list(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Ev_init, Im_init=Iv_init, H=H, M=M, t_vec=times, maxIterations=maxIterations)
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
#' @param Ev_init initial number of infected but not yet infectious mosquitoes.
#' @param Iv_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#'
#' @export

RM1_stochastic_sync <- function(a=0.3, p=0.9, mu=NULL, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Ev_init=0, Iv_init=0, H=100, M=100, times=0:100) {
	
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
	args <- list(a=a, mu=mu, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Em_init=Ev_init, Im_init=Iv_init, H=H, M=M, times=times)
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
#' @param Ev_init initial number of infected but not yet infectious mosquitoes.
#' @param Iv_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M_init initial mosquito population size (number of adult female mosquitoes).
#' @param times vector of times at which output should be returned.
#' @param Ktimes vector of times at which carrying capacity is defined.
#' @param Kvalues vector of carrying capacities that come into action at \code{Ktimes}.
#'
#' @export

RM2_deterministic <- function(a=0.3, p=0.9, mu=NULL, lambda=0.2, u=22, v=10, r=1/200, b=1, c=1, Eh_init=0, Ih_init=10, Ev_init=0, Iv_init=0, H=100, M_init=100, times=0:100, Ktimes=c(0,100), Kvalues=c(100,200)) {
	
	# calculate mu from p if not specified
	if (is.null(mu))
		mu <- -log(p)
	
	# solve ode
	mod <- RM2_deterministic_odin(a=a, mu=mu, lambda=lambda, u=u, v=v, r=r, b=b, c=c, Eh_init=Eh_init, Ih_init=Ih_init, Ev_init=Ev_init, Iv_init=Iv_init, H=H, M_init=M_init, Ktimes=Ktimes, Kvalues=Kvalues)
	output <- as.data.frame(mod$run(times))
	names(output)[1] <- 'time'
	
	return(output)
}
RM2_deterministic_odin <- odin::odin({
	
	# interpolation of carrying capacity
	K <- interpolate(Ktimes, Kvalues, "linear")
	
	# delay states
	Sh_du <- delay(Sh, u, 0)
	Iv_du <- delay(Iv, u, 0)
	Sv_dv <- delay(Sv, v, 0)
	Ih_dv <- delay(Ih, v, 0)

	# mosquito population rates of change
	deriv(M) <- lambda*M*(1-M/K) - mu*M
	
	# human flows between states
	EIR <- a*Iv/H
	FOI <- b*EIR
	
	Sh_to_Eh <- FOI*Sh
	Eh_to_Ih <- (a*b*Iv_du/H)*Sh_du # equivalent to delayed force of infection * delayed Sh
	Ih_to_Sh <- r*Ih
	
	# human rates of change
	deriv(Sh) <- -Sh_to_Eh + Ih_to_Sh
	deriv(Eh) <- Sh_to_Eh - Eh_to_Ih
	deriv(Ih) <- Eh_to_Ih - Ih_to_Sh
	
	# mosquito flows between states
	Sv_to_Ev <- a*c*Sv*Ih/H
	Ev_to_Iv <- a*c*Sv_dv*Ih_dv/H*exp(-mu*v)
	
	# mosquito rates of change
	deriv(Sv) <- -Sv_to_Ev - mu*Sv + lambda*M*(1-M/K)
	deriv(Ev) <- Sv_to_Ev - Ev_to_Iv - mu*Ev
	deriv(Iv) <- Ev_to_Iv - mu*Iv
	
	# initial conditions
	initial(M) <- M_init
	
	initial(Sh) <- H - Eh_init - Ih_init
	initial(Eh) <- Eh_init
	initial(Ih) <- Ih_init
	
	initial(Sv) <- M - Ev_init - Iv_init
	initial(Ev) <- Ev_init
	initial(Iv) <- Iv_init
	
	# parameters
	a <- user()
	mu <- user()
	lambda <- user()
	u <- user()
	v <- user()
	r <- user()
	b <- user()
	c <- user()
	Eh_init <- user()
	Ih_init <- user()
	Ev_init <- user()
	Iv_init <- user()
	H <- user()
	M_init <- user()
	Ktimes[] <- user()
	dim(Ktimes) <- user()
	Kvalues[] <- user()
	dim(Kvalues) <- user()
	
	# other outputs
	output(H) <- H
	output(K) <- K
	output(EIR) <- EIR
	output(Sh_to_Eh) <- Sh_to_Eh
	output(Eh_to_Ih) <- Eh_to_Ih
	output(Ih_to_Sh) <- Ih_to_Sh
	output(Sv_to_Ev) <- Sv_to_Ev
	output(Ev_to_Iv) <- Ev_to_Iv
})

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

