reglm <- function(formula = formula(data),
		random = NULL,
		family = gaussian(),
		data = sys.parent(),
		weights = NULL,
		dispersion = NULL,
		subset = NULL,
		na.action = na.fail,
		contrasts = NULL,
		start = NULL,
		mustart = NULL,
		betastart = NULL,
		sigmastart = NULL,
		control = reglm.control(...),
		ykeep = T,
		xkeep = F,
		ukeep = F,
		...)
{
#
#  Generalized linear models with random effects as in Schall, R (1991).
#  Estimation in glms with random effects. Biometrika, 78, 719-727.
#  Gordon Smyth
#  7 Oct 99.  Last revised 6 Nov 2004.
#
#  If no random effects, just use glm
	nrandom <- length(random)
	if(nrandom==0) {
		if(is.null(start) && !is.null(mustart)) start <- family()$link(mustart)
		return(glm(formula = formula,
			family = family,
			data = data,
			subset = subset,
			na.action = na.action,
			weights = weights,
			contrasts = contrasts,
			start = start,
			control = glm.control(...),
			model = model,
			y = ykeep,
			x = xkeep,
			...))
	}
#
#  Check that random effects are supplied as list
	if(!is.list(random)) stop("Random factors must be supplied as a list, e.g., list(A=a)")
#
#  Set up fixed effects model: 
#               y   response
#               X   design matrix
#               w   prior weights
#          offset   offset in linear predictor
#          family   response family
#           start  starting values for linear predictor (optional)
#         mustart   starting values for mu (optional)
#       betastart   starting values for coefficients (optional)
#
#	Save call for future reference
	scall <- match.call()
#
#	Evaluate the model frame
	mnames <- c("", "formula", "data", "weights", "subset", "na.action")
	cnames <- names(scall)
	cnames <- cnames[match(mnames,cnames,0)]
	mcall <- scall[cnames]
	mcall[[1]] <- as.name("model.frame")
	mframe <- eval(mcall, sys.parent())
#
#	Now extract the glm components
	y <- model.extract(mframe, response)
	N <- length(y)
	mterms <- attr(mframe, "terms")
	X <- model.matrix(mterms, mframe, contrasts)
	w <- model.extract(mframe, weights)
	if(is.null(w)) w <- rep(1,N)
	offset <- model.extract(mframe, offset)
	family <- as.family(family)
	p <- ncol(X)
	if(is.null(family$family["name"])) family$family["name"] <- "Unnamed"
#
#	Does this family have a known dispersion?
	fn <- family$family["name"]
	prior.dispersion <- dispersion
	if(is.null(dispersion) && (fn=="Binomial"||fn=="Poisson")) prior.dispersion <- 1
	if(is.null(dispersion)) dispersion <- 1
#
#  Set up random effects model:
#               U   design matrix 
#      sigmastart   starting values for variance components (optional)
#
#	Evaluate the model frame and extract components
	q <- rep(0,nrandom)
	for (i in 1:nrandom) {
		if(!is.factor(random[[i]])) random[[i]] <- factor(random[[i]])
		q[i] <- length(levels(random[[i]]))
	}
	qsum <- sum(q)
	U <- matrix(0,N,qsum)
	j0 <- 0
	for (i in 1:nrandom) {
		fval <- as.numeric(random[[i]])
		U[,j0+1:q[i]] <- fval%*%matrix(1,1,q[i])==matrix(1,N,1)%*%(1:q[i])
		j0 <- j0+q[i]
	}
#
#	Augmented design matrix (X2 is C)
	z2 <- rep(0,N+qsum)
	W2 <- rep(1,N+qsum)
	X2 <- matrix(0,N+qsum,p+qsum)
	qsumi <- 1:qsum
	X2[1:N,1:p] <- X
	X2[1:N,p+qsumi] <- U
	diag( X2[N+qsumi,p+qsumi] ) <- 1
#
#	Design matrix for the variance components
	Z <- rep(1:nrandom,q)%*%matrix(1,1,nrandom)==matrix(1,qsum,1)%*%matrix(1:nrandom,1,nrandom)
	mode(Z) <- "numeric"
#
#  Starting values for mean
	eval(family$initialize)
	if(!is.null(betastart)) {
		eta <- X %*% betastart
		mu <- family$inverse(eta+offset)
	}
	else if(!is.null(mustart)) {
		mu <- mustart
		eta <- family$link(mu)-offset
	}
	else if(!is.null(start)) {
		eta <- start
		mu <- family$inverse(eta+offset)
	}
	else eta <- family$link(mu)-offset
#
#	Starting value for variance components
	if(!is.null(sigmastart)) 
		sigma <- sigmastart
	else
		sigma <- rep(max(c(var(eta),1e-6)),nrandom)
#
#	Schall iteration
	epsilon <- control$epsilon
	maxit <- control$maxit
	trace <- control$trace
	iter <- 0
	repeat {
		iter <- iter+1
		if(iter>maxit) {
			warning("Max iterations exceed")
			break
		}
#
#		BLUP
		W <- eval(family$weight)
		z <- family$deriv(mu)*(y-mu)+eta
		z2[1:N] <- z
		W2[1:N] <- W/dispersion
		W2[N+qsumi] <- rep(1/sigma,q)
		fit <- lm.wfit(X2,z2,W2,method="qr")
		eta <- fit$fitted[1:N]
		mu <- family$inverse(eta+offset)
#
#		Update sigma
		rinv <- solve(fit$R, diag(p+qsum))
		varest <- rinv^2 %*% rep.int(1, p+qsum)
		bb <- t(Z) %*% coef(fit)[p+qsumi]^2
		v <- t(Z) %*% varest[p+qsumi]
		edf <- q - v/sigma
		sigmaold <- sigma
		sigma <- bb / edf
		if(is.null(prior.dispersion))
			dispersion <- sum(W*fit$residuals[1:N]^2) / (N-p-sum(edf))
		if(trace) cat("Iter",iter,sigma,dispersion,"\n")
		if( max(abs(sigma-sigmaold)/v/q) < epsilon) break
	}
	if(!is.null(prior.dispersion))
		dispersion <- sum(W*fit$residuals[1:N]^2) / (N-p-sum(edf))
#
#  Output for conditional model:
#  Exactly as for glm.object.  As for lm.object except that
#  linear.predictors and prior.weights are new components, and fitted.values
#  has a new definition.
	fit$family <- family$family
	fit$linear.predictors <- eta+offset
	fit$fitted.values <- mu
	if(!is.null(weights)) fit$prior.weights <- w
	fit$terms <- mterms
	fit$call <- scall
	fit$deviance <- family$deviance(mu,y,w)
	fit$null.deviance <- glm.null(X, y, w, offset, family)
	if(length(fit$null.deviance)>1) fit$null.deviance <- fit$null.deviance$null.deviance
	if(ykeep) fit$y <- y
	if(xkeep) fit$x <- X
	fit$formula <- as.vector(attr(mterms, "formula"))
#
#  Output for variance components:
	fit$sigma <- drop(sigma)
	fit$dispersion <- dispersion
	fit$prior.dispersion <- prior.dispersion
	fit$q <- q
	fit$edf <- drop(edf)
	fit$iter <- iter
	if(ukeep) fit$u <- U
	class(fit) <- c("glm","lm")
	fit
}

reglm.control <- function(epsilon = 1e-6, maxit = 50, trace = F, ...) {
#  Control iteration for reglm
#  As for glm.control but with different defaults
#  GKS  7 Oct 99
#
	if(epsilon <= 0) {
		warning("the value of epsilon supplied is zero or negative; the default value of 1e-7 was used instead")
		epsilon <- 1e-007}
	if(maxit < 1) {
		warning("the value of maxit supplied is zero or negative; the default value of 50 was used instead")
		maxit <- 50}
	list(epsilon = epsilon, maxit = maxit, trace = as.logical(trace)[1])
}
