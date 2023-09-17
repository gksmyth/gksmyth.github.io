remlscore <- function(y,X,Z,trace=F,tol=1e-5,maxit=40) {
#
#  Mean-variance fit by REML scoring
#  Fit normal(mu,phi) model to y with
#  mu=X%*%beta and log(phi)=Z%*%gam
#
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#  11 Sept 2000.  Last modified 28 Dec 2000.

n <- length(y)
p <- dim(X)[2]
q <- dim(Z)[2]
const <- n*log(2*pi)

# initial residuals from unweighted regression
fitm <- lm.fit(X,y,method="qr",singular.ok=T,qr=T)
Q <- qr.Q(fitm$qr)
h <- as.vector(Q^2 %*% array(1, c(p, 1)))
d <- fitm$residuals^2

# starting values
# use of weights guarantee that regression can be computed even if 1-h = 0
wd <- 1-h
zd <- log( d/(1-h) )+1.27
fitd <- lm.wfit(Z,zd,wd,method="qr",singular.ok=T)
gam <- ifelse(is.na(fitd$coef),0,fitd$coef)
g <- fitd$fitted.values
phi <- exp(g)
wm <- 1/phi
fitm <- lm.wfit(X,y,wm,method="qr",singular.ok=T,qr=T)
d <- fitm$residuals^2
dev <- sum(d/phi)+sum(log(phi))+const+2*log(prod(abs(diag(fitm$R))))

# reml scoring
iter <- 0
if(trace) cat("Iter =",iter,", Dev =",dev," Gamma",gam,"\n")
Q2 <- array(0,c(n,p*(p+1)/2))
repeat {
	iter <- iter+1

	# information matrix and leverages
	Q <- qr.qy(fitm$qr, diag(1, nrow = n, ncol = p))
	j0 <- 0
	for(k in 0:(p-1)) {
		Q2[ ,(j0+1):(j0+p-k)] <- Q[ ,1:(p-k)] * Q[ ,(k+1):p]
		j0 <- j0+p-k
	}
	Q2[ ,(p+1):(p*(p+1)/2)] <- sqrt(2) * Q2[ ,(p+1):(p*(p+1)/2)]
	h <- Q2[ ,1:p] %*% array(1,c(p,1))
	Q2Z <- t(Q2) %*% Z
	ZVZ <- ( t(Z) %*% vecmat(1-2*h,Z) + t(Q2Z) %*% Q2Z	)/2
	maxinfo <- max(diag(ZVZ))
	if(iter==1) {
		lambda <- abs(mean(diag(ZVZ)))/q
		I <- diag(q)
	}

	# score vector
	zd <- ( d - (1-h)*phi ) / phi
	dl <- crossprod(Z,zd)/2

	# Levenberg damping
	gamold <- gam
	devold <- dev
	lev <- 0
	repeat {
		lev <- lev+1

		# trial step
		dgam <- as.vector(solve.Hermitian( ZVZ + lambda*I, dl ))
		gam <- gamold + dgam
		phi <- as.vector(exp( Z %*% gam ))
		wm <- 1/phi
		fitm <- lm.wfit(X,y,wm,method="qr",singular.ok=T,qr=T)
		d <- fitm$residuals^2
		dev <- sum(d/phi)+sum(log(phi))+const+2*log(prod(abs(diag(fitm$R))))
		if(dev < devold - 1e-15) break

		# exit if too much damping
		if(lambda/maxinfo > 1e15) {
			gam <- gamold
			warning("Too much damping - convergence tolerance not achievable")
			break
		}

		# step not successful so increase damping
		lambda <- 2*lambda
		if(trace) cat("Damping increased to",lambda,"\n")
	}

	# iteration output
	if(trace) cat("Iter =",iter,", Dev =",dev," Gamma",gam,"\n")

	# keep exiting if too much damping
	if(lambda/maxinfo > 1e15) break

	# decrease damping if successful at first try
	if(lev==1) lambda <- lambda/10

	# test for convergence
	if( crossprod(dl,dgam) < tol ) break

	# test for iteration limit
	if(iter > maxit) {
		warning("reml: Max iterations exceeded")
		break
	}
}

# Nominal standard errors
se.gam <- sqrt(diag(solve.Hermitian(ZVZ)))

Rinv <- solve(fitm$R,diag(p))
se.beta <- sqrt(diag(Rinv %*% t(Rinv)))

list(beta=fitm$coef,se.beta=se.beta,gamma=gam,se.gam=se.gam,mu=fitm$fitted,phi=phi,dev=dev,h=h)
}
