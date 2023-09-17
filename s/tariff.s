tariff <- function(formula = formula(data),
		dformula = ~1,
		nclaims = NULL,
		exposure = stop("Exposure (number of units at risk) must be specified"),
		link.power = 0,
		dlink.power = 0,
		var.power = 1.5,
		data = sys.parent(),
		subset = NULL,
		contrasts = NULL,
		method = "ml",
		mustart = NULL,
		betastart = NULL,
		phistart = NULL,
		control = dglm.control(...),
		ykeep = T,
		xkeep = F,
		zkeep = F,
		...)
{
#
#	Model insurance claims data using Tweedie's compound Poisson model
#	Modified version of dglm
#	Gordon Smyth, Walter and Eliza Hall Institute.
#	30 Jan 99.  Revised 24 July 2001.
#
#	Check input arguments
#
	if(!is.null(nclaims))
		if(any(nclaims<0)) stop("Number of claims (nclaims) cannot be negative")
	if(any(exposure<0)) stop("Cannot have negative exposure")
	if(var.power <= 1 || var.power >= 2) stop("var.power must be strictly between 1 and 2")
#
#	Set up mean submodel: 
#               y   response - cost per unit as risk
#               X   design matrix
#          offset   offset in linear predictor
#         mustart   starting values for mu (optional)
#       betastart   starting values for coefficients (optional)
#
#	Save call for future reference
	scall <- match.call()
#
#	Evaluate the model frame
	scall$weights <- scall$exposure
	mnames <- c("", "formula", "data", "weights", "subset")
	cnames <- names(scall)
	cnames <- cnames[match(mnames,cnames,0)]
	mcall <- scall[cnames]
	mcall[[1]] <- as.name("model.frame")
	mframe <- eval(mcall, sys.parent())
	scall$weights <- NULL
#
#	Now extract the glm components
	y <- model.extract(mframe, response)
	n <- length(y)
	mterms <- attr(mframe, "terms")
	X <- model.matrix(mterms, mframe, contrasts)
	pw <- model.extract(mframe, weights)
	offset <- model.extract(mframe, offset)
	family <- tweedie(var.power=var.power, link.power=link.power)
	p <- var.power
	a <- (2-p)/(1-p)
#
#	Set up dispersion submodel:
#               Z   design matrix 
#         doffset   offset in linear predictor
#         dfamily   family for unit deviances
#        phistart   starting values for phi (optional)
#
#	Setup dformula but with y on left hand side
	mcall$formula <- formula
	mcall$formula[3] <- switch(match(length(scall$dformula),c(0,2,3)),
		1,scall$dformula[2],scall$dformula[3])
#
#	Evaluate the model frame and extract components
	mframe <- eval(mcall, sys.parent())
	dterms <- attr(mframe, "terms")
	Z <- model.matrix(dterms, mframe, contrasts)
	dpw <- rep(1,n)
	doffset <- model.extract(mframe, offset)
	dfamily <- tweedie(var.power=2, link.power=dlink.power)
#
#	Subset claims if necessary
	if(!is.null(subset) && !is.null(nclaims)) nclaims <- nclaims[subset]
#
#	Match method (ml or reml)
#
	if(!is.null(scall$method)) {
		name.method <- substitute(method)
		if(!is.character(name.method))
			name.method <- deparse(name.method)
		list.methods <- c("ml","reml","ML","REML","Ml","Reml")
		i.method <- pmatch(method,list.methods,nomatch=0)
		if(!i.method) stop("Method must be ml or reml")
		method <- switch(i.method,"ml","reml","ml","reml","ml","reml")
	}
	reml <- method=="reml"
#
#	Starting values for mu and phi
#
	y <- as.numeric(y)
	if(!is.null(betastart)) {
		eta <- X %*% betastart
		mu <- family$inverse(eta+offset)}
	else if(!is.null(mustart)) {
		mu <- mustart
		eta <- family$link(mu) - offset}
	else {
		mu <- y + (y==0) / 6
		eta <- family$link(mu) - offset}
	if(!is.null(phistart))
		phi <- phistart
	else
		phi <- rep(1,n)
#
#	First iteration of mean submodel
#	(necessary to get adjusted profile likelihood correct)
#
	gdot <- family$deriv(mu)
	zm <- eta + (y - mu) * gdot
	wm <- pw / mu^p / phi / gdot^2
	mfit <- lm.wfit(X, zm, wm, method="qr", singular.ok=T, qr=reml)
	eta <- mfit$fitted.values
	mu <- family$inverse(eta+offset)
#
#	Compute deviances.  Estimate constant phi if no phistart
#	(even if constant not in span of Z; necessary to match mean model iteration)
#
	if(is.null(nclaims)) {
		d <- family$deviance(mu, y, pw, residuals=T)^2
		if(is.null(phistart)) phi <- rep(mean(d),n)
		const <- glm.constant(y, family)}
	else {
		theta <- mu^(1-p) / (1-p)
		kappa <- mu^(2-p) / (2-p)
		t <- y * theta - kappa
		if(is.null(phistart)) phi <- rep( mean(-pw * t) / mean(nclaims) * (p-1), n)
		dpw <- 2 * pw * kappa / (p-1) / phi
		d <- -2 / dpw * ( pw * t + phi * nclaims / (p-1) ) + phi
		n1 <- (nclaims>0)
	}
	deta <- dfamily$link(phi)
#
#  Estimate model by alternate iterations
#
	epsilon <- control$epsilon
	maxit <- ceiling(control$maxit)
	trace <- control$trace
	iter <- 0
	m2loglik <- Inf
	repeat {
		iter <- iter+1
#
#		dispersion submodel
		gddot <- dfamily$deriv(phi)
		wd <- dpw / 2 / phi^2 / gddot^2
		if(reml) {
			h <- hat(mfit$qr)
			d <- d * dpw / (dpw-h)
			wd <- wd * pmax(dpw-h,0) / dpw
		}
		zd <- deta + (d - phi) * gddot
		dfit <- lm.wfit(Z, zd, wd, method="qr", singular.ok=T)
		deta <- dfit$fitted.values
		phi <- dfamily$inverse(deta+doffset)
#
#		mean submodel
		gdot <- family$deriv(mu)
		zm <- eta + (y - mu) * gdot
		wm <- pw / mu^p / phi / gdot^2
		mfit <- lm.wfit(X, zm, wm, method="qr", singular.ok=T, qr=reml)
		eta <- mfit$fitted.values
		mu <- family$inverse(eta+offset)
		if(is.null(nclaims))
			d <- family$deviance(mu, y, pw, residuals=T)^2
		else {
			theta <- mu^(1-p) / (1-p)
			kappa <- mu^(2-p) / (2-p)
			dpw <- 2 * pw * kappa / (p-1) / phi
			t <- y * theta - kappa
			d <- -2 / dpw * ( pw * t + phi * nclaims / (p-1) ) + phi
		}
#
#		overall likelihood
		m2loglikold <- m2loglik
		if(is.null(nclaims))
			m2loglik <- const + sum(log(phi/pw) + d/phi)
		else {
			m2loglik <- -2 * sum( pw*t/phi - lgamma(nclaims+1) )	- 2 * sum(
				-log(y[n1]) - lgamma(-nclaims[n1]*a) + nclaims[n1] * (
				log(pw[n1]/phi[n1])/(p-1) - a*log(y[n1]/(p-1)) - log(2-p) ) )}
		if(reml)
			m2loglik <- m2loglik + 2*log(abs(prod(diag(mfit$R))))

		if(trace)
			cat("DGLM iteration ", iter, ": -2*log-likelihood = ",
			format(round(m2loglik, 4)), " \n", sep = "")
		if ( abs(m2loglikold-m2loglik) < epsilon || iter == maxit )
			break
	}
#
#  Output for mean model:
#  Exactly as for glm.object.
#
	mfit$family <- family$family
	mfit$linear.predictors <- mfit$fitted.values+offset
	mfit$fitted.values <- mu
	mfit$prior.weights <- pw
	mfit$terms <- mterms
	mfit$call <- scall
	mfit$deviance <- family$deviance(mu, y, pw/phi)
	mfit$null.deviance <- glm.null(X, y, pw/phi, offset, family)
	if(length(mfit$null.deviance)>1) mfit$null.deviance <- mfit$null.deviance$null.deviance
	if(ykeep) mfit$y <- y
	if(xkeep) mfit$x <- X
	mfit$formula <- as.vector(attr(mterms, "formula"))
#
#  Output for dispersion model:
#  As for glm.object.  Nested in one output component.
#
#	In full likelihood case, set response to be unit estimators of phi
#	This is artificial because the actual response is bivariate
	if(!is.null(nclaims)) d <- -pw * t / pmax(nclaims,0.5) * (p-1)

	dfit$family <- dfamily$family
	dfit$linear.predictors <- dfit$fitted.values+doffset
	dfit$fitted.values <- phi
	if(!is.null(nclaims)) dfit$prior.weights <- dpw
	dfit$terms <- dterms
	scall$formula <- scall$dformula
	scall$dformula <- NULL
	scall$family <- call("tweedie",var.power=2,link.power=link.power)
	dfit$call <- scall
	dfit$deviance <- dfamily$deviance(phi, d, dpw/2)
	dfit$null.deviance <- glm.null(Z, d, dpw/2, doffset, dfamily)
	if(length(dfit$null.deviance)>1) dfit$null.deviance <- dfit$null.deviance$null.deviance
	if(ykeep) dfit$y <- d
	if(zkeep) dfit$z <- Z
	dfit$formula <- as.vector(attr(dterms, "formula"))
	dfit$iter <- iter
	class(dfit) <- c("glm","lm")

	out <- c(mfit, list(dispersion.fit = dfit, iter=iter, method=method, m2loglik=m2loglik))
	class(out) <- c("dglm","glm","lm")
	out
}
