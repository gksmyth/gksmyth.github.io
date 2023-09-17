pronyfreq <- function(y,nfreq=1,constant=T,tol=1e-5,maxit=10,trace=F,warnings=T) {
#	Estimate frequencies by constrained eigenvalue methods.
#	Includes constant term by default.

#	Smyth, G. K. (1999). Employing symmetry constaints for improved frequency
#	estimation by eigenanalysis methods. Submitted.
#	Kahn, M., Mackisack, M. S., Osborne, M. R., and Smyth, G. K. (1992). On
#	the consistency of Prony's method and related algorithms. J. Comput. Graph.
#	Statist. 1, 329-349.

#	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	2 July 99.  Last revised 10 Oct 99

# Check input
n <- length(y)
minn <- 4
maxf <- (n-1)/3
if(!constant) {
	minn <- minn-1
	maxf <- maxf+1/3
}
if(n < minn) stop(paste("Need at least",minn,"observations"))
if(nfreq > maxf) stop(paste("Can't estimate",nfreq,"frequencies"))

# Number of complex exponentials
p <- 2*nfreq+1
if(!constant) p <- p-1

# Form Q
Q <- matrix(0,p+1,nfreq+1)
sqrt05 <- 1/sqrt(2)
if(constant) {
	for (i in 1:(nfreq+1)) {
		Q[i,i] <- sqrt05
		Q[p+1-i+1,i] <- -sqrt05
	}
} else {
	for (i in 1:nfreq) Q[p+1-i+1,i] <- Q[i,i] <- sqrt05
	Q[nfreq+1,nfreq+1] <- 1
}

# form Y
i <- 0:(n-p-1)
Y <- matrix(0,n-p,p+1)
for (j in 1:(p+1)) Y[,j] <- y[j+i]
YQ <- Y %*% Q

# Constrained Pisarenko
B <- crossprod(YQ) / (n-p)
g <- eigen(B,symmetric=T)$vectors[,nfreq+1]
d <- Q %*% g

# Constrained ORA
nM <- as.integer(n-p)
mM <- as.integer(p+1)
M <- matrix(0,nM,mM)
mode(M) <- "double"
LYQ <- YQ
iter <- 0
repeat {
	iter <- iter+1
	if(iter > maxit) {
		if(warnings) warning("No convergence in maxit iterations")
		break
	}
	for (j in 1:(p+1)) M[,j] <- t(d[1:(p+2-j)]) %*% d[j:(p+1)]
	choleb <- .Fortran("choleb", nM, mM, LDL=M)
	D <- choleb$LDL[,1]
	if(any(D <= 0)) {
		if(warnings) warning("singularity encountered")
		break
	}
	D <- sqrt(D)
	for (j in 1:(nfreq+1)) {
		solblt <- .Fortran("solblt", nM, mM, LDL=choleb$LDL, LYQ=YQ[,j])
		LYQ[,j] <- solblt$LYQ / D
	}
	B <- crossprod(LYQ) / (n-p)
	gold <- g
	g <- eigen(B,symmetric=T)$vectors[,nfreq+1]
	d <- Q %*% g
	if(trace) cat("Iter",iter,round(g,6),"\n")
	if(max(abs((g-gold) / ifelse(g!=0,g,1e-5))) < tol) break
}
if(d[1]==0) return(NULL)
d <- as.vector(d)

# Extract known root = 1
if(constant) {
	d <- cumsum(d)[1:(2*nfreq+1)]
	d[(nfreq+2):(2*nfreq+1)] <- d[nfreq:1]
}

# Extract frequencies
freq <- Arg(polyroot(d))
freq <- unique( sort(freq)[(nfreq+1):(2*nfreq)] )

# Fit to data
r <- matrix(0:(n-1),n,1)
f <- matrix(freq,1,length(freq))
X <- cbind( cos( r%*%f ), sin( r%*%f ) )
if(constant) X <- cbind(1,X)
qrX <- qr(X)

list(prony=d,freq=freq,coef=qr.coef(qrX,y),fitted=qr.fitted(qrX,y),
     residuals=qr.resid(qrX,y),iter=iter)
}
