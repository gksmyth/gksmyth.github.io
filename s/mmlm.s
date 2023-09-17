library(Matrix)

mmlm <- function(X,y)
{
#	Multi-stage robust estimation for the linear model
#	(roughly equivalent to the S-Plus function of the same name)
#	X is the design matrix not including the intercept
#	GKS  18 Oct 01

	n <- dim(X)[1]
	p <- dim(X)[2]
	quan <- floor((n+p+1)/2)

#	Step 1. Approximate LTS regression using a genetic algorithm
	fit.lts <- ltsreg(X,y,singular.ok=T,quan=quan)

#	Step 2. LS regression on best subset of observations
	best <- ( abs(fit.lts$residuals) < quantile(abs(fit.lts$residuals),prob=(quan+0.5)/n) )
	Xbest <- cbind(1,X[best,])
	ybest <- y[best]
	fit.ls <- lm.fit(Xbest,ybest)

#	Step 3. M-estimation of scale
	scale <- mscale(fit.ls$residuals)

#	Step 4. M-estimation of mean
	X <- cbind(1,X)
	linpred <- function(X,b) X %*% b
	desmat <- function(X,b) X
	fit.mm <- mmnl(X,y,b=fit.ls$coef,fun=linpred,grad=desmat,scale=scale)
	fit.mm
}
