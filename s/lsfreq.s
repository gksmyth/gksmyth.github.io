lsfreq <- function(y,x=NULL,freq,constant=F,trace=F)
{
#	Frequency estimation by separable least squares
#
#	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	12 Jun 99.  Last revised 15 Oct 99.
#
	n <- length(y)
	nfreq <- length(freq)
	if(is.null(x)) x <- 0:(n-1)
	X1 <- matrix(1,n,constant+2*nfreq)
	xf <- x %*% matrix(freq,1,nfreq)
	cosi <- constant + 1:nfreq
	sini <- cosi + nfreq
	X1[,cosi] <- cos(xf)
	X1[,sini] <- sin(xf)
	qr1 <- qr(X1)
	z <- qr.resid(qr1,y)
	rho <- crossprod(z)
	if (rho/max(abs(y)) < 1e-15) return(freq=freq,coef=qr.coef(qr1,y),fitted=qr.coef(qr1,y),residuals=z)
#
	iter <- 0
	repeat {
		iter <- iter+1
		coef <- qr.coef(qr1,y)
		X2 <- vecmat(x,(matvec(cos(xf),coef[sini])-matvec(sin(xf),coef[cosi])))
		X2.1 <- qr.resid(qr1,X2)
		d1 <- crossprod(X2.1,z)
		d2 <- crossprod(X2.1)
		maxd2 <- max(d2)
		if(iter==1) {
			lambda <- sqrt(mean(d2^2))/1000
			I <- diag(nfreq)
		}
		# Levenberg damping
		freqold <- freq
		old <- rho
		lev <- 0
		repeat {
			lev <- lev+1
			dfreq <- solve.Hermitian( d2 + lambda*I, d1 )
			freq <- freqold + dfreq
			xf <- x %*% matrix(freq,1,nfreq)
			X1[,cosi] <- cos(xf)
			X1[,sini] <- sin(xf)
			qr1 <- qr(X1)
			z <- qr.resid(qr1,y)
			rho <- crossprod(z)
			if(rho < old - 1e-15 || rho==0) break
			if(lambda/maxd2 > 1e15) {
				freq <- freqold
				warning("Levenberg tolerance not achievable")
				break
			}
			lambda <- 2*lambda
			if(trace) cat("Lambda",lambda,"\n")
		}
		if(trace) cat("Iter",iter,"Freq",round(freq,5),"SS",round(rho,8),"\n")
		if(lambda/maxd2 > 1e15) break
		if(lev==1) lambda <- lambda/10
		# Test for convergence
		if( crossprod(d1,dfreq) < 1e-7 ) break
		if(iter > 50) {
			warning("Max iterations exceeded")
			break
		}
	}
	list(freq=as.vector(freq),coef=qr.coef(qr1,y),fitted=qr.fitted(qr1,y),residuals=z)
}
