
# -----------------------------------
header = function(name="Name") {

s <- paste("
# ",name,".R

# Author: Bob Verity
# Date: ",Sys.Date(),"

# Purpose:
# (this is an example header)

# ------------------------------------------------------------------

",sep="")
cat(s)

}

# -----------------------------------
myFunctions = function() {
    
s <- "
    #### FUNCTION LIST
    
    #### Misc
    # loadPackage
    # zeropad
    # logdescriptive
    # perlin_noise
    # colVars
    # rowVars
    # rateRatio
    # grad
    # cubeSpline_segment
    # cubeSpline
    # harmonicMean
    # first
    # last
    # latexTable
    # changeNames
    # reportConsole
    # dot
    # newLine
    # monthDays
    # is.int
    # removeCols
    # matrix.smooth
    # minSpanTree
    # fast
    # %+=% (increment by some value)
    # vec2mat
    # bin2D
    
    #### Probability
    # rdirichlet
    # rt_scaled
    # rMultiMix
    # normal_loglike
    # dBDI
    # dloggamma
    # rSTR1
    # dSTR1
    # rCRP
    # rlomax
    # dlomax
    # rDPM
    
    #### Combinatorics
    # partitions
    # StirlingFirst
    # StirlingSecond
    # increment_restrictedGrowth
    # increment_vector
    # allsamps
    # convertRadix
    
    #### Plotting
    # plotOn
    # plotOff
    # animateOn
    # animateOff
    # MCMCPlot
    # MCMCPlot2
    # MCMCBiPlot
    # densityPlot
    # errorBars
    # pieCharts
    # win
    # coordText
    # multiPanel
    # image_fix
    # filledContour2
    
    #### Colour Palettes
    # colPlot
    # transHex
    # bobRainbow
    # bobRedBlue
    # bobRedBlue2
    # smoothCols
"
cat(s)
    
}

# -----------------------------------
## plotOn()
plotOn = function(active=FALSE,fileroot='',filestem1='',filestem2='',fileindex='',type='pdf',width=6,height=6,res=300) {
#### start optional save-to-file function
	
	if (active) {
		if (type=='pdf') {
			pdf(file=paste(fileroot,filestem1,filestem2,fileindex,".pdf",sep=""), bg="white", width=width, height=height)
		}
		if (type=='png') {
			png(file=paste(fileroot,filestem1,filestem2,fileindex,".png",sep=""), bg="white", width=width*500, height=height*500,res=res)
		}
	}
}

# -----------------------------------
## plotOff()
plotOff = function(active=FALSE) {
#### end optional save-to-file function
	
	if (active) {
		dev.off()
	}
}

# -----------------------------------
## zeropad()
zeropad = function(number,padding) {
#### add leading zeros to number. Can handle negative numbers.
	
	preamble <- ''
	if (substr(number,1,1)=='-') {
		preamble <- '-'
		number <- substr(number,2,nchar(x))
	}
	output <- paste(paste(rep(0,padding-nchar(number)),collapse=""),number,sep="")
	output <- paste(preamble,output,sep='')
	return(output)
}

# -----------------------------------
## rdirichlet()
rdirichlet = function(alpha_vec) {
#### draw from Dirichlet distribution (non-symmetric)
	
	Y = rgamma(length(alpha_vec),shape=alpha_vec,scale=1)
	output = Y/sum(Y)
	return(output)
}

# -----------------------------------
## rt_scaled()
rt_scaled = function(n,df,ncp,scale) {
#### draw from scaled t distribution
	
	sigma = 1/rgamma(n,shape=df/2,rate=df*scale^2/2)
	output = rnorm(n,mean=ncp,sd=sqrt(sigma))
	return(output)
}

# -----------------------------------
## dt_scaled()
dt_scaled = function(x,df,ncp,scale,log=FALSE) {
#### density of scaled t distribution
	
	output = lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2)-((df+1)/2)*log(1+1/df*((x-ncp)/scale)^2)
	if (log==FALSE) output = exp(output)
	return(output)
}

# -----------------------------------
## normal_loglike()
normal_loglike = function(x,sd,priormean,priorsd) {
#### log-probability of data x integrated over normal likelihood with normal prior
	
	n = length(x)
	xbar = mean(x)
	S = sum((x-xbar)^2)
	v = n*priorsd^2/sd^2+1
	output = -n/2*log(2*pi*sd^2) - 0.5*log(v) - 1/(2*sd^2)*(S+n*(xbar-priormean)^2/v)
	if (output=="NaN") output=0
	return(output)
}

# -----------------------------------
## logDescriptive()
logDescriptive = function(Y) {
### Calculates descriptive stats on extremely large or extremely small values. Input vector of values in log-space. Output sample mean and variance, also in log-space.
	
	n = length(Y)
	log_ybar = log(mean(exp(Y-min(Y))))+min(Y)
	log_samplevar = 2*min(Y) + log(n/(n-1)*mean(exp(2*Y-2*min(Y))-2*exp(Y+log_ybar-2*min(Y))+exp(2*log_ybar-2*min(Y))))
	return(c(log_ybar,log_samplevar))
}

# -----------------------------------
## perlin_noise()
perlin_noise = function(outRows=100,outCols=100,levelsX=10,levelsY=10) {
#### generates perlin noise of any scale in a matrix of any size (code copied from internet)
#### outRows = rows in output matrix
#### outCols = columnds in output matrix
#### levelsX = bumpyness in x dimension
#### levelsY = bumpyness in y dimension
	
	# convert to more convenient names
	M = outRows
	N = outCols
	n = levelsX
	m = levelsY
	
  # For each point on this n*m grid, choose a unit 1 vector
  vector_field <- apply(
    array( rnorm( 2 * n * m ), dim = c(2,n,m) ),
    2:3,
    function(u) u / sqrt(sum(u^2))
  )
  f <- function(x,y) {
    # Find the grid cell in which the point (x,y) is
    i <- floor(x)
    j <- floor(y)
    stopifnot( i >= 1 || j >= 1 || i < n || j < m )
    # The 4 vectors, from the vector field, at the vertices of the square
    v1 <- vector_field[,i,j]
    v2 <- vector_field[,i+1,j]
    v3 <- vector_field[,i,j+1]
    v4 <- vector_field[,i+1,j+1]
    # Vectors from the point to the vertices
    u1 <- c(x,y) - c(i,j)
    u2 <- c(x,y) - c(i+1,j)
    u3 <- c(x,y) - c(i,j+1)
    u4 <- c(x,y) - c(i+1,j+1)
    # Scalar products
    a1 <- sum( v1 * u1 )
    a2 <- sum( v2 * u2 )
    a3 <- sum( v3 * u3 )
    a4 <- sum( v4 * u4 )
    # Weighted average of the scalar products
    s <- function(p) 3 * p^2 - 2 * p^3
    p <- s( x - i )
    q <- s( y - j )
    b1 <- (1-p)*a1 + p*a2
    b2 <- (1-p)*a3 + p*a4
    (1-q) * b1 + q * b2
  }
  xs <- seq(from = 1, to = n, length = N+1)[-(N+1)]
  ys <- seq(from = 1, to = m, length = M+1)[-(M+1)]
  outer( xs, ys, Vectorize(f) )
}

# -----------------------------------
## dBDI()
dBDI = function(x,lambda,mu,alpha,log=FALSE) {
#### probability mass function for the birth-death-immigration model
	
	output1 = (alpha/lambda)*log(1-lambda/mu)
	output2 = x*log(lambda) + lgamma(x+alpha/lambda)-lgamma(alpha/lambda)-lgamma(x+1)-x*log(mu) + (alpha/lambda)*log(1-lambda/mu)
	output = (x==0)*output1 + (x>0)*output2
	if (!log)
		output = exp(output)
	return(output)
}

# -----------------------------------
## dlgamma()
dlgamma = function(x,alpha,beta,log=FALSE) {
#### density function for log(X), where X is gamma(shape=alpha,rate=beta)
	
	output = alpha*x+alpha*log(beta)-lgamma(alpha)-beta*exp(x)
	if (log==FALSE)
		output = exp(output)
	return(output)
}

# -----------------------------------
## partitions()
partitions = function(n) {
#### return all unique partitions of the numbers 1:n. These numbers can then be used to index a list if neccessary.
	
	if (n==1) {
		output = 1
	}
	if (n>1) {
		output = rbind(c(1,1),c(1,2))
	}
	if (n>2) {
		for (i in 3:n) {
			maxVals = mapply(max,split(output,f=1:nrow(output)))+1
			newSeq = mapply(seq,to=maxVals)
			output = cbind(output[rep(1:nrow(output),times=maxVals),],unlist(newSeq))
		}
	}
	return(output)
}

# -----------------------------------
## StirlingFirst()
StirlingFirst = function(n) {
#### function for computing the Stirling numbers of the first kind (output in log space)
	
	if (n<3) {
		output = rep(0,n)
	} else {
		Stirlingvec = c(1,1)
		K = 0
		for (i in 2:(n-1)) {
			tempmat = rbind(c(0,Stirlingvec),i*c(Stirlingvec,0))
			Stirlingvec = colSums(tempmat)
			K = K + log(max(Stirlingvec))
			Stirlingvec = Stirlingvec/max(Stirlingvec)
		}
		output = K + log(Stirlingvec)
	}
	return(output)
}

# -----------------------------------
## StirlingSecond()
StirlingSecond = function(n) {
#### function for computing the Stirling numbers of the second kind
	
	if (n==1) {
		output = 1
	} else {
		Stirlingvec = 1
		for (i in 1:(n-1)) {
			Stirlingvec = c(0,Stirlingvec) + (1:(i+1))*c(Stirlingvec,0)
		}
		output = Stirlingvec
	}
	return(output)
}

# -----------------------------------
## increment_restrictedGrowth()
increment_restrictedGrowth = function(group,K) {
#### find next group that satisfies restricted growth function (max value K). Allows user to eplore all unique partitions by iteratively updating the grouping. Returns 0 when no more groups.
	
	target = length(group)
	group[target] = group[target] + 1
	while (group[target]>(max(group[1:(target-1)])+1) | group[target]>K) {
		group[target] = 1
		target = target - 1
		group[target] = group[target] + 1
	}
	if (target==1) {
		return(0)
	} else {
		return(group)
	}
}

# -----------------------------------
## allsamps()
allsamps = function(targetlist) {
#### output all possible sequences that are possible by sampling (without replacement) from a given list
	
	lengths = mapply(length,targetlist)
	outputmat = matrix(0,nrow=prod(lengths),ncol=length(lengths))
	for (i in 1:length(lengths)) {
		val1 = prod(lengths[1:i][-i])
		val2 = prod(lengths[i:length(lengths)][-1])
		outputmat[,i] = rep(rep(targetlist[[i]],each=val2),times=val1)
		}
	return(outputmat)
	}

# -----------------------------------
## convertRadix()
convertRadix = function(x,base=2,fixlength=NA) {
#### convert decimal number system to any other radix (base)
	
	floornumber = x
	output = NULL
	while (floornumber!=0) {
		output = c(floornumber%%base,output)
		floornumber = floor(floornumber/base)
	}
	if (!is.na(fixlength)) {
		output = c(rep(0,fixlength-length(output)),output)
	}
	return(output)
}

# -----------------------------------
## dSTR1()
dSTR1 = function(k,n,theta,log=FALSE) {
#### probability mass function of STR1 (Stirling type 1) distribution. Defined for k in 1:n
	
	output = lgamma(theta)-lgamma(n+theta)+StirlingFirst(n)[k]+k*log(theta)
	if (!log) output = exp(output)
	return(output)
}

# -----------------------------------
## rSTR1()
rSTR1 = function(n,size,theta) {
### draw from STR1 (Stirling type 1) distribution
	
	drawmat = matrix(runif(size*n),nrow=size,ncol=n)
	output = colSums(drawmat<theta/(row(drawmat)-1+theta))
	return(output)
}

# -----------------------------------
## rSTR1()
rCRP = function(n,theta) {
#### draw group and group frequencies from a Chinese Restaurant Process
	
	if (n==1) return(list(group=1,groupFreqs=1))
	group = rep(1,n)
	maxGroup = 1
	groupFreqs = rep(0,n)
		groupFreqs[1] = 1
	groupFreqs2 = groupFreqs
		groupFreqs2[maxGroup+1] = theta
	for (i in 2:n) {
		newGroup = sample(n,1,prob=groupFreqs2)
		group[i] = newGroup
		if (newGroup>maxGroup) {
			maxGroup = newGroup
			groupFreqs2[maxGroup] = groupFreqs[maxGroup] = 1
			groupFreqs2[maxGroup+1] = theta
		} else {
			groupFreqs[newGroup] = groupFreqs[newGroup]+1
			groupFreqs2[newGroup] = groupFreqs2[newGroup]+1
		}
	}
	return(list(group=group,groupFreqs=groupFreqs))
}

# -----------------------------------
## colVars()
colVars = function(mat,sampleVar=TRUE,mu=0) {
#### calculate variance (sample or population) of columns of a matrix or data frame. If sampleVar==FALSE then calcualte relative to mean mu
	
	n = colSums(!is.na(mat))
	X = colSums(mat,na.rm=T)
	X2 = colSums(mat^2,na.rm=T)
	if (sampleVar) {
		output = 1/(n-1)*(X2-X^2/n)
	} else {
		output = 1/n*(X2-2*mu*X+n*mu^2)
	}
	return(output)
}

# -----------------------------------
## rowVars()
rowVars = function(mat,sampleVar=TRUE,mu=0) {
#### calculate variance (sample or population) of rows of a matrix or data frame. If sampleVar==FALSE then calcualte relative to mean mu
	
	n = rowSums(!is.na(mat))
	X = rowSums(mat,na.rm=T)
	X2 = rowSums(mat^2,na.rm=T)
	if (sampleVar) {
		output = 1/(n-1)*(X2-X^2/n)
	} else {
		output = 1/n*(X2-2*mu*X+n*mu^2)
	}
	return(output)
}

# -----------------------------------
## rlomax()
rlomax = function(n,alpha,lambda) {
#### draws from the lomax (Pareto type II) distribution with shape alpha and scale lambda
	
	r = rgamma(n,shape=alpha,rate=lambda)
	output = rexp(n,rate=r)
	return(output)
}

# -----------------------------------
## dlomax()
dlomax = function(n,alpha,lambda) {
#### density function for the lomax (Pareto type II) distribution with shape alpha and scale lambda
	
	output = log(alpha)+alpha*log(lambda)-(alpha+1)*log(x+lambda)
	if (log==FALSE) output = exp(output)
	return(output)
}

# -----------------------------------
## rateRatio()
rateRatio <- function(count1,time1,count2,time2,alpha=0.05) {
#### function for computing upper and lower confidence intervals on a ratio of rates
	
	point <- time2/time1*count1/count2
	LL <- time2/time1*count1/(count2+1)*qf(0.025,2*(count2+1),2*count1)
	UL <- time2/time1*(count1+1)/count2*qf(1-0.025,2*(count1+1),2*count2)
	if (count1==0 | time1==0 | count2==0 | time2==0) {
		point=LL=UL=NaN
	}
	return(list(point=point,LL=LL,UL=UL))
}

# -----------------------------------
## rMultiMix()
rMultiMix = function(n,K,psi,nu,alpha,tau,pMissing=0) {
#### draw from multivariate normal mixture model, with option for missing data
#### n = total number of draws
#### K = number of mixture components
#### psi = prior matrix for inverse Wishart prior on covariance
#### nu = degrees of freedom for inverse Wishart prior on covariance
#### alpha = Dirichlet prior parameter for mixture weights
#### tau = standard deviation of prior on components means (mean=0)
#### pMissing = proportion of data that is missing
	
	# load package mvtnorm for multivariate normal draws
	loadPackage('mvtnorm')
	
	# load package LaplacesDemon for Wishart normal draws
	loadPackage('LaplacesDemon')
	
	# generate missing data matrix (if possible, otherwise error)
	d = ncol(psi)
	pMissing = round(n*d*pMissing)/(n*d)
	if (pMissing<=0) {
		missingMat = NULL
	} else {
		if (pMissing>((d-1)/d)) {
			stop(paste("cannot generate ",pMissing*100,"% missing data without eliminating some observations entirely. Maximum missing data for this dimensionality is ",100*(d-1)/d,"%.",sep=""))
			missingMat = NULL
		} else {
			# draw number of missing elements in each observation
			missingRowSums = rmultinom(1,n*d*pMissing,prob=rep(1,n))
			while(max(missingRowSums)>(d-1)) {
				numOverLimit = sum(missingRowSums>(d-1))
				missingRowSums[missingRowSums>(d-1)] = missingRowSums[missingRowSums>(d-1)]-1
				numUnderLimit = sum(missingRowSums<(d-1))
				missingRowSums[missingRowSums<(d-1)] = missingRowSums[missingRowSums<(d-1)] + rmultinom(1,numOverLimit,prob=rep(1,numUnderLimit))
			}
			# create a random matrix that satisfies this constraint
			missingMat = matrix(0,n,d)
			for (j in 1:d) {
				probvec = missingRowSums/(d-j+1)
				missingMat[,j] = as.integer(runif(n)<probvec)
				missingRowSums = missingRowSums-missingMat[,j]
			}
		}
	}
	
	# draw from mixture model
	weights = rdirichlet(rep(alpha/K,K))
	groups = rep(1:K,times=rmultinom(1,n,prob=weights))
	draws = NULL
	sigma = list()
	mu = list()
	for (i in 1:K) {
		if (sum(groups==i)>0) {
			sigma[[i]] = rinvwishart(nu=nu,S=psi)
			mu[[i]] = rmvnorm(1,sigma=diag(tau^2,d))
			draws = rbind(draws,rmvnorm(sum(groups==i),mean=mu[[i]],sigma=sigma[[i]]))
		}
	}
	
	# return list of elements
	outlist = list(draws=draws,weights=weights,groups=groups,sigma=sigma,mu=mu,missingMat=missingMat)
	return(outlist)
}

# -----------------------------------
## grad()
grad = function(y,x=1:length(y)) {
#### numerically calculate gradient of values (x,y)
	
	delta_y = y[-1]-y[-length(x)]
	delta_x = x[-1]-x[-length(x)]
	grad = delta_y/delta_x
	grad_left = c(grad[1],grad)
	grad_right = c(grad,grad[length(grad)])
	grad_mean = (grad_left+grad_right)/2
	return(grad_mean)
}

# -----------------------------------
## cubeSpline_segment()
cubeSpline_segment = function(x,y,g) {
#### calculate cubic spline of the form ax^3+bx^2+cx+d between the points (x[1],y[1]) and (x[2],y[2]) with gradients g[1] and g[2]. Return constants a, b, c and d.
	
	# calculate coefficients
	a = (y[2]-y[1]-0.5*(g[1]+g[2])*(x[2]-x[1]))/(2*x[1]^3-3*x[1]^2*x[2]+x[2]^3-0.5*(3*x[2]^2-3*x[1]^2)*(x[2]-x[1]))
	b = (g[2]-g[1]-a*(3*x[2]^2-3*x[1]^2))/(2*x[2]-2*x[1])
	c = g[1]-3*a*x[1]^2-2*b*x[1]
	d = y[2]-a*x[2]^3-b*x[2]^2-g[1]*x[2]+3*a*x[1]^2*x[2]+2*b*x[1]*x[2]
	
	# return results
	return(list('a'=a,'b'=b,'c'=c,'d'=d))
}

# -----------------------------------
## cubeSpline()
cubicSpline = function(x,y,g,interpoints=10) {
#### calculate cubic spline passing through the points (x,y) with known gradient g at these points. The parameter 'interpoints' describes the number of points in each segment.
	
	# create smooth vectors
	n = length(x)
	xvec = seq(x[1],x[n],l=(n-1)*interpoints)
	yvec = rep(0,(n-1)*interpoints)
	
	# interpolate each segment using cubic spline
	for (i in 2:n) {
		whichPoints = ((i-2)*interpoints+1):((i-1)*interpoints)
		z = cubeSpline_segment(x[(i-1):i],y[(i-1):i],g[(i-1):i])
		yvec[whichPoints] = z$a*xvec[whichPoints]^3+z$b*xvec[whichPoints]^2+z$c*xvec[whichPoints]+z$d
	}
	
	return(list('x'=xvec,'y'=yvec))
}

# -----------------------------------
## harmonicMean()
harmonicMean = function(z) {
#### calculates harmonic mean of values z (already in log space), accounting for underflow. Returns value in log space.

	A = min(z)
	return(A-log(sum(exp(-(z-A))))+log(length(z)))	
}

# -----------------------------------
## MCMCPlot()
MCMCPlot = function(x, col=grey(0.7), pch=4, cex=0.4, xmin=1, xmax=length(x), ymin=min(x,na.rm=T), ymax=max(x,na.rm=T), xlab='iteration', main='', ylab=main) {
#### produces simple plot useful for visualising MCMC chains

	plot(x,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=main,pch=pch,cex=cex,col=col)
}

# -----------------------------------
## MCMCBiPlot()
MCMCBiPlot = function(x, y, col=grey(0.7), pch=4, cex=0.4, xmin=min(x,na.rm=T), xmax=max(x,na.rm=T), ymin=min(y,na.rm=T), ymax=max(y,na.rm=T), xlab='', ylab='', main='') {
#### produces bivariate plot useful for visualising correlation between two MCMC chains

	plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=main,pch=pch,cex=cex,col=col)
}

# -----------------------------------
## densityPlot()
densityPlot = function(x, xlab='', ylab='density', main='', boxOn=TRUE, meanOn=TRUE, medianOn=TRUE, xmin=1.5*min(x,na.rm=T)-0.5*max(x,na.rm=T), xmax=1.5*max(x,na.rm=T)-0.5*min(x,na.rm=T), ymin=0, ymax='auto') {
#### produce kernel density plot of data x, with mean, median and 95% limits shown
	
	# produce kernel density
	fx = density(x,from=xmin,to=xmax)
	if (ymax=='auto')
		ymax <- max(fx$y)
	plot(fx,ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=main)
	polygon(c(fx$x,rev(fx$x)), c(fx$y,rep(0,length(fx$y))), col=grey(0.7))
	
	# add quantiles to density plot
	quantiles = quantile(x,prob=c(0.025,0.975))
	fx_left_x = fx$x[fx$x<quantiles[1]]
	fx_left_y = fx$y[fx$x<quantiles[1]]
	polygon(c(fx_left_x,rev(fx_left_x)), c(fx_left_y,rep(0,length(fx_left_y))), col=colors()[373])
	fx_right_x = fx$x[fx$x>quantiles[2]]
	fx_right_y = fx$y[fx$x>quantiles[2]]
	polygon(c(fx_right_x,rev(fx_right_x)), c(fx_right_y,rep(0,length(fx_right_y))), col=colors()[373])
	
	# add lines for median and mean
	if (medianOn) {
		best_median = which.min(abs(fx$x-median(x)))
		lines(rep(fx$x[best_median],2),c(0,fx$y[best_median]))
	}
	if (meanOn) {
		best_mean = which.min(abs(fx$x-mean(x)))
		lines(rep(fx$x[best_mean],2),c(0,fx$y[best_mean]),lty=2)
	}
	
	# add legend
	if (boxOn) {
		legend('topright', legend=paste(c('   mean=','median=',' q0.025=',' q0.975='), signif(c(mean(x),median(x), quantiles[1],quantiles[2])), sep=''), lty=c(2,1,NA,NA), seg.len=1, bg='#FFFFFF80', box.col='#00000050')
	}
}

# -----------------------------------
## bobRainbow()
bobRainbow = function() {
#### 12 rainbow colours
	
	# use function rgb(red/255, green/255, blue/255) to convert RGB to hex
	#ee2c2c=RGB(238,44,44)
	#64b4ff=RGB(100,180,255)
	#7ce42c=RGB(124,228,44)
	#ffdc00=RGB(255,220,0)
	#aa46ff=RGB(170,70,255)
	#ff7f00=RGB(255,127,0)
	#ff50aa=RGB(255,80,170)
	#2c962c=RGB(44,150,44)
	#1f78c8=RGB(31,120,200)
	#be6e32=RGB(190,110,50)
	#6b6b6b=RGB(107,107,107)
	#c8b4ff=RGB(200,180,255)
	
	colours = c('#ee2c2c','#64b4ff','#7ce42c','#ffdc00','#aa46ff','#ff7f00','#ff50aa','#2c962c','#1f78c8','#be6e32','#6b6b6b','#c8b4ff')
	return(colours)
}

# -----------------------------------
## bobRedBlue()
bobRedBlue = function() {
#### 6 colours going from red to blue
	
	# use function rgb(red/255, green/255, blue/255) to convert RGB to hex
	#FB2D0A=RGB(251,45,10)
	#FED00B=RGB(254,208,11)
	#62D500=RGB(98,213,0)
	#1F7C1A=RGB(31,124,26)
	#1B79FE=RGB(27,121,254)
	#00006D=RGB(0,0,109)
	
	colours = c('#FB2D0A','#FED00B','#62D500','#1F7C1A','#1B79FE','#00006D')
	return(colours)
}

# -----------------------------------
## bobRedBlue2()
bobRedBlue2 = function() {
#### 11 colours going from red to blue
	
	# use function rgb(red/255, green/255, blue/255) to convert RGB to hex
	#9E0142=RGB(158,1,66)
	#D53E4F=RGB(213,62,79)
	#F46D43=RGB(244,109,67)
	#FDAE61=RGB(253,174,97)
	#FEE08B=RGB(254,224,139)
	#FFFFBF=RGB(255,255,191)
	#E6F598=RGB(230,245,152)
	#ABDDA4=RGB(171,221,164)
	#66C2A5=RGB(102,194,165)
	#3288BD=RGB(50,136,189)
	#5E4FA2=RGB(94,79,162)
	
	colours = c('#9E0142', '#D53E4F', '#F46D43', '#FDAE61', '#FEE08B', '#FFFFBF', '#E6F598', '#ABDDA4', '#66C2A5', '#3288BD', '#5E4FA2')
	return(colours)
}

# -----------------------------------
## loadPackage()
loadPackage = function(string) {
#### check for package and install/load if needed
#### for example, if string="devtools" the function will run;
#if (!"devtools"%in%rownames(installed.packages())) {
#	install.packages('devtools')
#}
#require("devtools")
	
	# create command string
	s1 = paste('if (!\'',string,'\'%in%rownames(installed.packages())) {',sep='')
	s2 = paste('install.packages(\'',string,'\')',sep='')
	s3 = '}'
	s4 = paste('require(',string,')',sep='')
	
	# evaluate string
	fullString = paste(s1,s2,s3,s4,sep='\n')
	eval(parse(text=fullString))
}

# -----------------------------------
## colPlot()
colPlot = function(colVec, size=8) {
#### produce simple plot to visualise a color palette
	
	n = length(colVec)
	plot(1:n, rep(0,n), col=colVec, pch=20, cex=size, axes=FALSE, xlab='', ylab='')
	text(1:n,rep(0,n),labels=1:n)
}

# -----------------------------------
## transHex()
transHex = function(c, alpha) {
#### take vector of colours c in hexadecimal form and build in transparancy alpha
	
	z = t(col2rgb(c))/255
	output = mapply(function(x) {rgb(x[1],x[2],x[3],alpha=alpha)}, split(z,f=row(z)))
	return(output)
}


# -----------------------------------
## first()
first = function(x) {
#### return first value in a vector

	x[1]
}

# -----------------------------------
## last()
last = function(x) {
#### return last value in a vector

	x[length(x)]
}

# -----------------------------------
## latexTable()
latexTable = function(x, digits=3, file='~/Desktop/tester.txt') {
#### convert x to LaTeX format with set decimal places and write to file

	outVec = rep(NA,nrow(x))
	s = paste("%.",digits,"f",sep='')
	for (i in 1:nrow(x)) {
		outVec[i] = paste(paste(sprintf(s,x[i,]),collapse=' & '),'\\\\')
	}
	write.table(outVec, file=file, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# -----------------------------------
## smoothCols()
smoothCols <- function(x, xmin=min(x,na.rm=T), xmax=max(x,na.rm=T), res=1e3, rawCols='DEFAULT') {
#### read in continuous values between xmin and xmax and return colours associated with these values, taken from a smoothly varying scale. Default colours are taken from tim.colors().
	
	# load RColorBrewer if needed
	loadPackage('RColorBrewer')
	
	# get x between 0 and 1
	x = (x-xmin)/(xmax-xmin)
	
	# define colours if needed
	if (rawCols[1]=='DEFAULT')
		rawCols = c('#00008F', '#00009F', '#0000AF', '#0000BF', '#0000CF', '#0000DF', '#0000EF', '#0000FF', '#0010FF', '#0020FF', '#0030FF', '#0040FF', '#0050FF', '#0060FF', '#0070FF', '#0080FF', '#008FFF', '#009FFF', '#00AFFF', '#00BFFF', '#00CFFF', '#00DFFF', '#00EFFF', '#00FFFF', '#10FFEF', '#20FFDF', '#30FFCF', '#40FFBF', '#50FFAF', '#60FF9F', '#70FF8F', '#80FF80', '#8FFF70', '#9FFF60', '#AFFF50', '#BFFF40', '#CFFF30', '#DFFF20', '#EFFF10', '#FFFF00', '#FFEF00', '#FFDF00', '#FFCF00', '#FFBF00', '#FFAF00', '#FF9F00', '#FF8F00', '#FF8000', '#FF7000', '#FF6000', '#FF5000', '#FF4000', '#FF3000', '#FF2000', '#FF1000', '#FF0000', '#EF0000', '#DF0000', '#CF0000', '#BF0000', '#AF0000', '#9F0000', '#8F0000', '#800000')
		
	# make smooth colours
	myPal = colorRampPalette(rawCols)
	cols = myPal(res+1)[floor(x*res)+1]
	return(cols)
}

# -----------------------------------
## changeNames()
changeNames <- function(df,oldNames,newNames) {
#### change data frame column names

	names(df)[names(df)%in%oldNames] <- newNames
	return(df)
}

# -----------------------------------
## reportConsole()
reportConsole <- function(i,imax,istep=1) {
#### use within loop to print current iteration to console
	
	if (i%%istep==0) {
		s <- paste("iteration ",i," of ",imax,", (",round(i/imax*100),"%)\n",sep="")
		cat(s)
		flush.console()
	}
}

# -----------------------------------
## dot()
dot <- function(n=1) {
#### prints one or more dots to screen. Follow by newLine() to add carriage return

	cat(paste(rep('.',n),collapse=''))
	flush.console()
}

# -----------------------------------
## newLine()
newLine <- function(n=1) {
#### prints one or more carriage returns

	cat(paste(rep('\n',n),collapse=''))
	flush.console()
}

# -----------------------------------
## monthDays()
monthDays <- function(year) {
#### return days in each month of chosen year

	date1 <- paste(year,"-01-01",sep="")
	date2 <- paste(year+1,"-01-01",sep="")
	output <- as.numeric(diff(seq(as.Date(date1),as.Date(date2),by="month")))
	return(output)
}

# -----------------------------------
## is.int()
is.int <- function(x) {
#### check that character string can be converted to integer without issue

	output = FALSE
	if (!is.na(suppressWarnings(as.numeric(x))))
		output = (as.numeric(x)%%1==0)
	return(output)
}

# -----------------------------------
## errorBars()
errorBars <- function(y1, y2, x=1:length(y1), se=FALSE, width=1, col=1, lty=1, lwd=1) {
#### produce error bars at given positions. Use se=TRUE to input mean and standard error, otherwise input raw upper and lower limits.

	# calculate limits based on method
	if (se) {
		LL <- y1-2*y2
		UL <- y1+2*y2
	} else {
		LL <- y1
		UL <- y2
	}

	segments(x,LL,x,UL,lty=lty,lwd=lwd,col=col)
	segments(x-width/2,LL,x+width/2,LL,lty=lty,lwd=lwd,col=col)
	segments(x-width/2,UL,x+width/2,UL,lty=lty,lwd=lwd,col=col)
}

# -----------------------------------
## increment_vector()
increment_vector <- function(x,binMax=c(1,4,3,2)) {
#### increment vector x under binMax constraint. Return NULL at final increment.

	if (length(x)!=length(binMax))
		stop("x and binMax of different lengths")
	if (any(x>binMax))
		stop("some x greater than binMax")
	if (all(x==binMax))
		return()
	for (i in length(x):1) {
		x[i] <- x[i]+1
		if (x[i]<=binMax[i]) {
			break
		} else {
			x[i] <- 1
		}
	}
	return(x)
}

# -----------------------------------
## removeCols()
removeCols <- function(df,nameVec) {
#### remove named columns from data frame

	output <- subset(df,select=setdiff(names(df),nameVec))
	return(output)
}

# -----------------------------------
## matrix.smooth()
matrix.smooth <- function(M,reps=1) {
#### interpolate rows and columns of a matrix alternately using spline. reps is the number of times to double the matrix size.

	for (i in 1:reps) {
		M2 <- cbind(M,M)
		for (i in 1:nrow(M)) {
			M2[i,] <- spline(M[i,],n=2*length(M[i,]))$y
		}
		M <- rbind(M2,M2)
		for (i in 1:ncol(M2)) {
			M[,i] <- spline(M2[,i],n=2*length(M2[,i]))$y
		}
	}
	return(M)
}

# -----------------------------------
## minSpanTree()
minSpanTree <- function(data) {
#### use Prim's algorithm to calculate minimum spanning tree. The 'data' argument must be a matrix (not a data frame) with observations in rows and dimensions in columns. Distances are calculated as Euclidean over any number of dimensions. Returns multiple objects; 1) a data frame of nodePairs giving all pairwise links that make up the tree in the order that they were added, 2) an edgeList giving the nodes attached to each node, 3) edgeNum giving the number of edges associated with each node.
	
	# extract basic properties of data
	n <- nrow(data)
	dims <- ncol(data)
	
	# initialise objects for storing results
	nodePairs <- data.frame(node1=rep(0,n-1),node2=0)
	edgeList <- replicate(n,NULL)
	edgeNum <- rep(0,n)
	
	# initialise algorithm by measuring all distances relative to first point
	point1 <- data[1,]
	data <- data[-1,]
	d_running <- 0
	for (j in 1:dims) {
		d_running <- d_running + (data[,j]-point1[j])^2
	}
	
	# given a point A and a cluster B, the distance between A and B is equal to the *minimum* distance between A and any point in B. The vector bestNode stores *which* node in B the point A links to
	bestNode <- rep(1,n-1)
	# during the algorithm the data object is trimmed, and so the row index no longer stores the unique ID of the data point. The vector pos fixes this by maintaining a record of unique IDs.
	pos <- 2:n
	
	# run Prims algorithm
	for (i in 1:(n-1)) {
		# node2 is the node with the minimum distance to the current tree. node1 is the node within the tree that node2 links to
		nodeIndex <- which.min(d_running)
		node1 <- bestNode[nodeIndex]
		node2 <- pos[nodeIndex]
		
		# store nodes and edges
		nodePairs$node1[i] <- node1
		nodePairs$node2[i] <- node2
		edgeList[[node1]] <- c(edgeList[[node1]], node2)
		edgeList[[node2]] <- c(edgeList[[node2]], node1)
		edgeNum[node1] <- edgeNum[node1]+1
		edgeNum[node2] <- edgeNum[node2]+1
		
		# calculate new distance of all points from node2
		d_new <- 0
		for (j in 1:dims) {
			d_new <- d_new + (data[,j]-data[nodeIndex,j])^2
		}
		# recalculate the *minimum* distance between the current tree and all nodes, and if the new distance is shorter then reflect that in bestNode
		bestNode[d_running>d_new] <- node2
		d_running <- (d_running<=d_new)*d_running + (d_running>d_new)*d_new
		
		# drop new node from data set and all other objects
		data <- data[-nodeIndex,,drop=FALSE]
		d_running <- d_running[-nodeIndex]
		pos <- pos[-nodeIndex]
		bestNode <- bestNode[-nodeIndex]
	}
	# return multiple objects
	return(list('nodePairs'=nodePairs, 'edgeList'=edgeList, 'edgeNum'=edgeNum))
}

# -----------------------------------
## rDPM()
rDPM <- function(n, sigma=1, priorMean_x=0, priorMean_y=0, priorSD=10, alpha=1) {
#### draw from a two-dimensional Dirichlet process mixture model. Return multiple outputs, including true group allocation
    
    # draw grouping
	group = rep(1,n)
	freqs = 1
	for (i in 2:n) {
		group[i] = sample(length(freqs)+1,1,prob=c(freqs,alpha))
		if (group[i]>length(freqs))
        freqs = c(freqs,0)
		freqs[group[i]] = freqs[group[i]] + 1
	}
	group = sort(group)
	
	# draw means and data
	mu_x = rnorm(length(freqs),priorMean_x,priorSD)
	mu_y = rnorm(length(freqs),priorMean_y,priorSD)
	x = rnorm(n,mu_x[group],sigma)
	y = rnorm(n,mu_y[group],sigma)
	
	# return results
	return(list(x=x, y=y, group=group, mu_x=mu_x, mu_y=mu_y))
}

# -----------------------------------
## pieCharts()
pieCharts <- function(x, y, proportions, radius=0.2, lwd=1, border_col=1, seg_col=1:max(mapply(length,proportions)), smoothness=51, x_stretch=1) {
#### adds pie charts to existing plot. proportions need not sum to one.

	theta <- seq(0,2*pi,l=smoothness)
	for (i in 1:length(x)) {
		segs <- length(proportions[[i]])
		z <- c(0, cumsum(proportions[[i]]/sum(proportions[[i]])))
		if (segs==1) {
			polygon(x[i]+x_stretch*radius*cos(theta), y[i]+radius*sin(theta), lwd=lwd, border=border_col, col=seg_col[1])
		} else {
			for (j in 2:(segs+1)) {
				theta_sub <- c(2*pi*z[j-1], theta[theta>(2*pi*z[j-1]) & theta<(2*pi*z[j])], 2*pi*z[j])
				polygon(c(x[i],x[i]+x_stretch*radius*sin(theta_sub),x[i]), c(y[i],y[i]+radius*cos(theta_sub),y[i]), lwd=lwd, border=border_col, col=seg_col[j-1])
			}
		}
	}
}

# -----------------------------------
## win()
win = function(rows=1,cols=1) {
#### shortcut for running par(mfrow=c(rows,cols)), which takes slightly longer to type!
	
	par(mfrow=c(rows,cols))
}

# -----------------------------------
## animateOn()
animateOn <- function(active=FALSE,fileName,saveFrames=FALSE,width=480,height=480,units="px",pointsize=12,quality=75,bg="white",res=NA) {
#### first of two functions needed to create mpg from static images. Run this function, then create a bunch of plots, then run animateOff() to close. animateOff() should have the same shared arguments as animateOn().
#### Note that this function requires ImageMagick to be installed to work, and has only been tested on Windows. If running for the first time then you will likely need to add ImageMagick to the path using something like the following:
#### Sys.setenv(PATH=paste(Sys.getenv("PATH"),"C:\\Program Files\\ImageMagick-6.9.3-Q16",sep=";"))
	
	if (!active)
		return()
	filePath <- paste(fileName,"%04d.jpeg",sep="")
	if (!saveFrames)
		filePath <- paste(tempdir(),"/",filePath,sep="")
	jpeg(filePath,width=width,height=height,units=units,pointsize=pointsize,quality=quality,bg=bg,res=res)
}

# -----------------------------------
## animateOff()
animateOff <- function(active=FALSE,fileName,saveFrames=FALSE,frameRate=25) {
#### second of two functions needed to create mpg from static images. This function should be run after running animateOn() and creating a series of plots. animateOff() should have the same shared arguments as animateOn().

	if (!active)
		return()
	dev.off()
	filePath <- paste(fileName,"%04d.jpeg",sep="")
	if (!saveFrames)
		filePath <- paste(tempdir(),"/",filePath,sep="")
	myCommand <- paste("ffmpeg -y -r ",frameRate," -i ",filePath," ",fileName,".mp4",sep="")
	shell(myCommand,intern=TRUE)
}

# -----------------------------------
## coordText()
coordText <- function(x=0.05,y=0.95,text='A)',cex=1.5) {
#### defines a square region from 0-1 in both x and y and inserts text at designated location
#### useful for e.g. adding panel labels to figures
	
	pars <- par(new=T,xpd=T,mar=c(0,0,0,0))
	plot(0, type='n', axes=F, xlab=NA, ylab=NA, xlim=c(0,1), ylim=c(0,1))
	text(x,y,text,cex=cex)
	par(pars)
}

# -----------------------------------
## MCMCPlot2()
MCMCPlot2 <- function(x, h_cells=100, v_cells=100, xmin=1, xmax=length(x), ymin=min(x,na.rm=T), ymax=max(x,na.rm=T), xlab='iteration', main='', ylab=main) {
#### similar to MCMCPlot, but for very large numbers of samples where ordinary plotting would not be feasible. Instead, bins points into cells before plotting.

	# make x,y data frame
	df <- data.frame(x=1:length(x),y=x)

	# make vectors of breaks and matrix z for storing final results
	v_breaks <- seq(ymin, ymax, l=v_cells+1)
	h_breaks <- seq(xmin, xmax, l=h_cells+1)
	v_mids <- (v_breaks[-1]+v_breaks[-length(v_breaks)])/2
	h_mids <- (h_breaks[-1]+h_breaks[-length(h_breaks)])/2
	z <- matrix(0,v_cells,h_cells)
	
	# split data based on vertical breaks
	df_split <- split(df,f=cut(df$y,v_breaks))
	
	# for each element, count frequency classes based on horizontal breaks
	for (i in 1:length(df_split)) {
		df_split2 <- split(df_split[[i]],f=cut(df_split[[i]]$x,h_breaks))
		z[i,] <- mapply(nrow,df_split2)
	}
	
	
	# plot final matrix
	palette <- colorRampPalette(colors()[c(1,131,107)])
	image(h_mids, v_mids, t(z), col=palette(100), xlab=xlab, ylab=ylab, main=main)
}

# -----------------------------------
## fast()
fast <- function(fileName, header='auto') {
#### shortcut function for reading in data quickly using the data.table package
	
	loadPackage('data.table')
	return(as.data.frame(fread(fileName, header=header)))
}

# -----------------------------------
## multiPanel()
multiPanel <- function() {
#### this function is not really a proper function (it does not return anything). Rather, it is a place to store this useful bit of code that can adapted as needed. Copy this code over, and then drop whatever individual plots are needed into the inner loop.

# setup multi-panel plotting parameters
plot_rows <- 2
plot_cols <- 2
outer_margin <- c(4,4,2,1)
tickLength <- -0.5
x_range <- c(-5,5)
y_range <- c(0,300)
x_axis_at <- seq(-4,4,2)
y_axis_at <- seq(0,250,50)
sub_main <- paste('submain',1:(plot_rows*plot_cols))
x_lab <- 'x axis'
y_lab <- 'y axis'
x_cex <- 0.8
y_cex <- 0.8

# create multi-panel plot
par_store <- par(mfrow=c(plot_rows, plot_cols), mar=c(0,0,0,0), oma=outer_margin, tcl=tickLength)
index <- 0
for (i in 1:plot_rows) {
	for (j in 1:plot_cols) {
		index <- index+1
		
		# put individual plots here
		hist(rnorm(1e3), xlim=x_range, ylim=y_range, ann=FALSE, axes=FALSE)
		title(sub_main[index],line=-1)
		
		# add axes and box
		if (i==plot_rows)
			axis(1, at=x_axis_at)
		if (j==1)
			axis(2, at=y_axis_at)
		box()
	}
}
mtext(x_lab, side=1, outer=TRUE, cex=x_cex, line=2.2)
mtext(y_lab, side=2, outer=TRUE, cex=y_cex, line=2.2)
par(par_store)

}

# -----------------------------------
## %+=%
# infix operator for incrementing a variable by a user-defined amount
'%+=%' <- function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

# -----------------------------------
## image_fix()
image_fix <- function(z,x=1:ncol(z),y=1:nrow(z), ...) {
#### produces an image() plot, but flips input matrix so that plot has same orientation as raw object
	
	image(x,y,t(z[nrow(z):1,]), ...)
}

# -----------------------------------
## filledContour2()
filledContour2 <- function(z,x=NULL,y=NULL,l=11,col=bobRedBlue2(),fix=TRUE, ...) {
#### produces nice alternative version of filled contour
	
	# set x and y based on matrix dimensions
	x_auto <- is.null(x)
	y_auto <- is.null(y)
	if (x_auto)
		x <- 1:nrow(z)
	if (y_auto)
		y <- 1:ncol(z)
	
	# if fix orientation
	if (fix) {
		z <- t(z[nrow(z):1,])
		if (x_auto)
			x <- 1:ncol(z)
		if (y_auto)
			y <- 1:nrow(z)
	}
	
	# produce plot
	myLevels <- seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE), l=l+1)
	myCols <- smoothCols(1:l,rawCols=col)
	image(x,y,z, col=myCols)
	contour(x,y,z, levels=myLevels, drawlabels=FALSE, add=TRUE)
}

# -----------------------------------
## vec2mat()
vec2mat <- function(x,y,dim) {
#### produce matrix from two vectors. dim=1 returns the x-matrix, dim=2 the y-matrix
	
	if (dim==1) {
		output <- matrix(rep(x,each=length(y)),length(y))
	} else {
		output <- matrix(rep(y,length(x)),length(y))
	}
	return(output)
}

# -----------------------------------
## bin2D()
bin2D <- function(x,y,x_breaks,y_breaks) {
#### read in x and y data, along with vectors of breaks. Output 2D bin counts and vectors of midpoints.

	# bin data in both x and y
	freq <- as.data.frame(table(findInterval(x,x_breaks),findInterval(y,y_breaks)))
	freq[,1] <- as.numeric(as.character(freq[,1]))
	freq[,2] <- as.numeric(as.character(freq[,2]))
	
	# get rid of bins outside of range
	freq <- freq[freq[,1]>0 & freq[,1]<length(x_breaks) & freq[,2]>0 & freq[,2]<length(y_breaks),]
	
	# fill in 2D matrix
	freq2D <- matrix(0,length(y_breaks)-1,length(x_breaks)-1)
	freq2D[cbind(freq[,1],freq[,2])] <- freq[,3]
	
	# calculate midpoints
	x_mids <- (x_breaks[-1]+x_breaks[-length(x_breaks)])/2
	y_mids <- (y_breaks[-1]+y_breaks[-length(y_breaks)])/2
	
	# output all as list
	output <- list(x_mids= x_mids,y_mids= y_mids,z=freq2D)
}