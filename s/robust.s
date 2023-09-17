#	Distribution license: The software code in this file may be
#	freely copied, translated or redistributed providing that the
#	orginal source is acknowledged and the following publication
#	is cited in redistributed documentation:

#	Smyth, G. K., and Hawkins, D. M. (2000).
#	Robust frequency estimation using elemental sets.
#	Journal of Computational and Graphical Statistics 9, 196-214.


mmnl <- function(X,y,b,fun,grad,scale=NULL,trace=F)
{
#	Nonlinear MM-estimator with 50% breakdown and 95% efficiency
#	Based on methods described by Yohai (Annals of Statistics, 1987) and Stromberg (JASA, 1993)
#
#	Gordon Smyth
#	Created 3 Jun 99.  Last revised 15 Oct 1999.
#
	z <- (y-fun(X,b))
	if(all(z==0)) {
		if(is.null(scale)) scale <- 0
		return(b=b,fitted=y,residuals=rep(0,length(y)),scale=scale,criterion=0)
	}
	if(is.null(scale)) scale <- mscale(z)
	z <- z/0.9014/scale
	iter <- 0
	repeat {
		iter <- iter+1
		G <- grad(X,b)/0.9014/scale
		w <- wt.hampel(z,a=1.5,b=3.5,c=8)
		wG <- vecmat(sqrt(w),G)
		d1 <- crossprod(G,psi.hampel(z))
		d2 <- crossprod(wG)
		maxd2 <- max(d2)
		if(iter==1) {
			if(all(w==0)) stop("Specified scale too small")
			rho <- sum(rho.hampel(z))
			lambda <- sqrt(mean(d2^2))/1000
			I <- diag(dim(G)[2])
		}
		# Levenberg damping
		bold <- b
		old <- rho
		lev <- 0
		repeat {
			lev <- lev+1
			db <- solve.Hermitian( d2 + lambda*I, d1 )
			b <- bold + db
			z <- (y-fun(X,b))/0.9014/scale
			rho <- sum(rho.hampel(z))
			if(rho < old - 1e-15 || rho==0) break
			if(lambda/maxd2 > 1e15) {
				b <- bold
				warning("Levenberg tolerance not achievable")
				break
			}
			lambda <- 2*lambda
			if(trace) cat("Lambda",lambda,"\n")
		}
		if(trace) cat("Iteration",iter,"\nb",b,"\nObjective Function",rho,"\n")
		if(lambda/maxd2 > 1e15) break
		if(lev==1) lambda <- lambda/10
		# Test for convergence
		if( crossprod(d1,db) < 1e-8 ) break
		if(iter > 40) {
			warning("mmnl: Max iterations exceeded")
			break
		}
	}
	mu <- fun(X,b)
	list(b=b,fitted=mu,residuals=y-mu,scale=scale,criterion=rho)
}

mscale <- function(u)
{
#	Scale M-estimator with 50% breakdown
#	Yohai (1987) Annals, Stromberg (1993) JASA.
#
#	GKS  2 June 99
#
	if(mean(u==0) >= 0.5) return(0)
	U <- abs(u)
	s <- median(U)/0.6744898
	iter <- 0
	repeat {
		iter <- iter+1
		z <- u/0.212/s
		d1 <- mean(rho.hampel(z))-3.75
		d2 <- mean(z*psi.hampel(z))
		s <- s*(1+d1/d2)
		if(iter > 50) {
			cat("mscale: Max iterations exceeded")
			break
		}
		if(abs(d1/d2) < 1e-14) break
	}
	s	
}

rho.hampel <- function(u, a = 1.5, b = 3.5, c = 8)
{
#	Integral of Hampel's redescending psi function (Hampel, Ronchetti,
#	Rousseeuw and Stahel, 1986, Robust Statistics, Wiley, page 150).
#	Default values are as in Stromberg (1993) JASA.
#
#	GKS  31 May 99
#
	U <- abs(u)
	A <- (U <= a)	#increasing
	B <- (U > a) & (U <= b)	#flat
	C <- (U > b) & (U <= c)	#descending
	D <- (U > c)	# zero
	rho <- U
	rho[A] <- (U[A] * U[A])/2
	rho[B] <- a * (U[B] - a/2)
	rho[C] <- a * (b - a/2) + a * (U[C] - b) * (1 - (U[C] - b)/(c - b)/2)
	rho[D] <- (a * (b - a + c))/2
	rho
}

psi.hampel <- function(u, a = 1.5, b = 3.5, c = 8)
{
#	Hampel's redescending psi function (Hampel, Ronchetti,
#	Rousseeuw and Stahel, 1986, Robust Statistics, Wiley, page 150).
#	Default values are as in Stromberg (1993) JASA.
#
#	GKS  2 June 99
#
	U <- abs(u)
	B <- (U > a) & (U <= b)	#flat
	C <- (U > b) & (U <= c)	#descending
	D <- (U > c)	# zero
	psi <- u
	psi[B] <- sign(u[B]) * a
	psi[C] <- sign(u[C]) * a * (c - U[C])/(c - b)
	psi[D] <- 0
	psi
}

mmfreq <- function(y,x=NULL,freq,coef=NULL,constant=F,scale=NULL,trace=F)
{
#	MM-frequency estimator with 50% breakdown and 95% efficiency
#	Yohai (1987) Annals, Stromberg (1993) JASA.
#
#	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	2 Jul 99.  Last revised 15 Oct 99
#
	if(is.null(x)) x <- 0:(length(y)-1)
	if(any(is.na(y))) {
		x <- x[!is.na(y)]
		y <- na.omit(y)
	}
	nfreq <- length(freq)
	if(!is.null(coef)) if(length(coef) != constant+2*nfreq) stop("Need 2 coefficients for each frequency")
	if(length(y) < constant+3*nfreq) stop("Need at least 3 observations for each frequency")
	fun <- function(x,b) {
		nfreq <- length(b) %/% 3
		constant <- length(b) %% 3
		f <- matrix(b[1:nfreq],1,nfreq)
		cosi <- constant + (nfreq+1):(2*nfreq)
		sini <- cosi + nfreq
  		Xb <- cos(x%*%f) %*% b[cosi] + sin(x%*%f) %*% b[sini]
		if(constant) Xb <- Xb + b[nfreq+1]
		as.vector(Xb)
	}
	grad <- function(x,b) {
		nfreq <- length(b) %/% 3
		constant <- length(b) %% 3
		f <- matrix(b[1:nfreq],1,nfreq)
		cosi <- constant + (nfreq+1):(2*nfreq)
		sini <- cosi + nfreq
		G <- x*( matvec(cos(x%*%f),b[sini]) - matvec(sin(x%*%f),b[cosi]) )
		if(constant)
			return( cbind( G, 1, cos(x%*%f), sin(x%*%f) ) )
		else
			return( cbind( G, cos(x%*%f), sin(x%*%f) ) )
	}
	if(is.null(coef)) {
		f <- matrix(freq,1,nfreq)
		G <- cbind( cos(x%*%f), sin(x%*%f) )
		coef <- ltsreg(G,y,intercept=constant)$coefficients
	}
	out <- mmnl(x,y,c(freq,coef),fun,grad,scale=scale,trace=trace)
	out$freq <- Arg(exp(out$b[1:nfreq]*1i))
	out$coef <- out$b[(nfreq+1):(3*nfreq+constant)]
	out$b <- NULL
	out
}

robfreq <- function(y,nfreq=1,s=NULL,breakdown=0.5,trace=F,maxregs=6000) {
#	Robust frequency estimation.
#	Should have nearly 50% breakdown and 95% efficiency when no outliers are
#	present.  Uses an elemental set LTS estimator as starting value for an
#	MM-estimator.
#
#	y     - input time series
#	nfreq - number of frequencies to find
#	s     - error standard deviation
#	trace - write working elemental estimates to screen?
#	maxregs - maximum number of criterion evaluations allowed

#	Author: Gordon Smyth
#	Created 15 Oct 1999

n <- length(y)
x <- 0:(n-1)
if(breakdown > 0.5) breakdown <- 0.5
if(breakdown < 0) breakdown <- 0
trim <- floor(n*(1-breakdown))+floor(3*nfreq/2)
cosi <- 1:nfreq
sini <- cosi + nfreq

# Make sure data is not trivial
if(all(y==0)) return(list(freq=0,coef=c(0,0),fitted=y,residuals=y,scale=0,criterion=0))
if(var(y)==0) {
	m <- mean(y)
	return(list(freq=0,coef=c(m,0),fitted=rep(m,n),residuals=rep(0,n),scale=0,criterion=0))
}

# Prony ORA estimator
if(trace) cat(" PronyFreq")
ora <- pronyfreq(y,nfreq,constant=F,maxit=5,warnings=F)
e <- ora$residuals^2
tsmin <- sum( sort(e)[1:trim] )
omegats <- ora$freq
alphats <- ora$coef
ets <- e
if(trace) {
	cat(" (")
	cat(round(omegats,4),sep=",")
	cat(")")
}

# Check for exact fit
if(max(abs(ora$residual))/max(abs(y)) < 1e-11) return(c(ora,list(scale=0,criterion=0)))

# If low breakdown wanted, return immediately
if(trim >= n) return(ora)

# Number of steps or spacings in elemental set
steps <- 3*nfreq-1

# Stromberg's recommendation for number of element sets
stromberg <- ceiling( log(0.001)/log(1-breakdown^(3*nfreq)) )
if(stromberg < n) stromberg <- n
if(stromberg > maxregs) stromberg <- maxregs

# Find maxs corresponding to Stromberg
nelemsets <- s <- 0
repeat {
	s <- s+1
	new <- n - steps*s
	if(new < 1) {
 		smax <- s-1
		prob <- rep(1,smax)
		break
	}
	if(s>1) {
		nelemsetsthin <- nelemsets + new * (1+(s-2)*nfreq) / (1+(s-1)*nfreq)
		if(nelemsetsthin >= stromberg) {
			smax <- s
			prob <- rep(1,smax)
			prob[smax] <- (1+(s-2)*nfreq) / (1+(s-1)*nfreq)
			break
		}
	}
	nelemsets <- nelemsets + new
	if(nelemsets >= stromberg) {
		smax <- s
		prob <- rep(1,smax)
		break
	}
}

# Make sure widest elemental sets at least at least 7% of data range
# If smax is increased, thin out to required number of elemental sets
if(new > 1) {
	mins <- ceiling(n/14/steps)
	if(smax < mins) {
		scut <- smax
		smax <- mins
		s <- 1:smax
		prob <- pmin( 1, (1+(scut-1)*nfreq) / (1+(s-1)*nfreq) )
	}
	s <- 1:smax
	nsets <- n-steps*s
	prob <- prob*stromberg/sum(nsets*prob)
}

# Make sure no more than maxregs regressions, including harmonics
# If necessary, thin out to maximum number of regressions
s <- 1:smax
nregs <- (n-steps*s) * (1+(s-1)*nfreq)
if(sum(nregs*prob) > maxregs) {
	scut <- smax
	repeat {
		scut <- scut-1
		prob <- pmin( 1, (1+(scut-1)*nfreq) / (1+(s-1)*nfreq) )
		if(sum(nregs*prob) <= maxregs) {
			prob <- pmin( 1, (1+scut*nfreq) / (1+(s-1)*nfreq) )
			prob <- prob/sum(nregs*prob)*maxregs
			break
		}
		if(scut == 1) {
			prob <- prob/sum(nregs*prob)*maxregs
			break
		}
	}
}

# Try elemental estimators
Xebase <- matrix(0,3*nfreq,2*nfreq)
X <- Xbase <- matrix(0,n,2*nfreq)
i <- 0:steps

# Step through all spacings
if(trace) cat(" Elemental")
for (spacing in s) {
	if(trace) cat(" ",spacing,sep="")

#	Sets within each spacing
	nsets <- n-(3*nfreq-1)*spacing
	x1sample <- sample(nsets,round(prob[spacing]*nsets))
	for (x1 in x1sample) {
		xe <- x1+spacing*i
		ye <- y[xe]
		eomega <- elemfreq(ye,nfreq=nfreq)
		if ( !is.null(eomega) ) {

#			Set omega to base frequency
			eomegabase <- eomega/spacing

#			Compute criteria at base frequency
			eo <- matrix(eomegabase,1,nfreq)
			xoe <- (xe-1)%*%eo
			Xebase[,cosi] <- cos(xoe)
			Xebase[,sini] <- sin(xoe)
			xo <- x%*%eo
			Xbase[,cosi] <- cos(xo)
			Xbase[,sini] <- sin(xo)
			qrXe <- qr(Xebase)
			alpha <- qr.coef(qrXe,ye)
			mu <- Xbase %*% alpha
			e <- (y-mu)^2
			tsbase <- sum( sort(e)[1:trim] )
#			cat(spacing,x1,"0 0",tsbase,"\n")
			eomegats <- eomegabase
			ealphats <- alpha
			eets <- e
			tsminh <- tsbase
	
			if (spacing>1) {
#				Try to find a better combination of harmonics for this elemental set
				Xemin <- Xebase
				Xmin <- Xbase
				for (f in 1:nfreq) {
					Xe <- Xemin
					X <- Xmin
#					Search harmonics of frequency f
					for (h in 1:(spacing-1)) {
						eomegah <- acos(cos(eomegabase[f]+2*pi*h/spacing))
						xoe <- (xe-1)*eomegah
						xo <- x*eomegah
						Xe[,f] <- cos(xoe)
						Xe[,f+nfreq] <- sin(xoe)
						X[,f] <- cos(xo)
						X[,f+nfreq] <- sin(xo)
						qrXe <- qr(Xe)
						alpha <- qr.coef(qrXe,ye)
						mu <- X %*% alpha
						e <- (y-mu)^2
						ts <- sum( sort(e)[1:trim] )
#						cat(spacing,x1,f,h,ts,"\n")
						if (ts < tsminh) {
							tsminh <- ts;
							eomegats[f] <- eomegah;
							ealphats <- alpha;
							eets <- e
							Xemin <- Xe;
							Xmin <- X
						}
					}
				}
			}

#			Now have best harmonic. Test again global minimum.
			if (tsminh < tsmin) {
				tsmin <- tsminh
				omegats <- eomegats
				alphats <- ealphats
				ets <- eets
				if(trace) {
					cat(" (")
					cat(round(omegats,4),sep=",")
					cat(")")
				}
			}
		}
	}
}

# Terminate if all frequencies are zero
if(all(omegats==0))
	return(list(freq=omegats,coef=alphats,fitted=y-ets,residuals=ets))

# If some but not all frequencies are zero, replace with constant term
constant <- any(omegats==0)
if(constant) omegats <- omegats[omegats != 0]

# Least squares on best half of data
good <- sort(sort.list(e)[1:trim])
ls.ts <- lsfreq(y=y[good],x=x[good],freq=omegats,constant=constant)
if(trace) {
	cat(" LS (")
	cat(round(ls.ts$freq,4),sep=",")
	cat(")")
}

# Compute scale M estimate
s <- mscale(ls.ts$residuals)
if(trace) cat(" s=",round(s,3),sep="")

# Terminate if scale is zero
if(s/max(ls.ts$coef) < 1e-12) return(c(ls.ts,list(scale=s)))

# Location M estimate
mm <- mmfreq(y,freq=ls.ts$freq,coef=ls.ts$coef,scale=s,constant=constant)
mm.ora <- mmfreq(y,freq=ora$freq,coef=ora$coef,scale=s)
if(mm.ora$criterion < mm.ora$criterion) mm <- mm.ora
if(trace) {
	cat(" MM (")
	cat(round(mm$freq,4),sep=",")
	cat(")\n")
}
mm
}

elemfreq <- function(y,nfreq=NULL) {
#	Fit frequencies to elemental set
#	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	12 Jul 99.  Last revised 15 Oct 99.

if(is.null(nfreq)) nfreq <- length(y) %/% 3
if (nfreq==1) {
	if(y[2]==0) return(NULL)
	d2 <- -(y[1]+y[3])/y[2]
	if (abs(d2)>2 ) return(NULL) else return(acos(-d2/2))
}

i <- 0:(nfreq-1)
i1 <- 1
i2 <- 2*nfreq+1
b <- y[i1+i] + y[i2+i]
B <- matrix(0,nfreq,nfreq)
for (j in 1:(nfreq-1)) {
	i1 <- i1+1
	i2 <- i2-1
	B[,j] <- y[i1+i] + y[i2+i]
}
B[,nfreq] <- y[nfreq+1+i]
qrB <- qr(B)
if (qrB$rank < nfreq) return(NULL)
d <- -qr.coef(qrB,b)
f <- log(polyroot( c(1,d,d[(nfreq-1):1],1) ))
if (any(abs(Re(f))>1e-8))
	return(NULL)
else
	return(sort(Im(f))[(nfreq+1):(2*nfreq)])
}
