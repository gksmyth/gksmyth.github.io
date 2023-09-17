dpoisgam <- function(y, mu = stop("No mean specified"), phi = 1, p = 1.5)
{
#  Density function for the Tweedie family with 1 <= p <= 2.
#  Gordon Smyth, 7 Dec 98
#
#  Check input arguments are admissable
#
	if(any(mu < 0)) stop("Mean of Tweedie family must be non-negative")
	if(any(phi < 0))
		stop("Dispersion of Tweedie family must be non-negative")
	if(length(p) > 1)
		stop("Index of Tweedie family must be scalar")
	if(p < 1 | p > 2) stop("Index of Tweedie family must be between 1 and 2")
#
#  Make sure arguments are the same length
#
	ny <- length(y)
	nmu <- length(mu)
	nphi <- length(phi)
	n <- max(c(ny, nmu, nphi))
	if(ny < n) y <- rep(y, length = n)
	if(nmu < n) mu <- rep(mu, length = n)
	if(nphi < n) phi <- rep(phi, length = n)
#
#  The Tweedie distribution can be expressed as a Poisson(lambda) sum of
#  gamma(alpha,beta)
#
	lambda <- mu^(2 - p)/(2 - p)/phi
	alpha <- (2 - p)/(p - 1)
	beta <- phi * (p - 1) * mu^(p - 1)
#
#  Special cases
#
	if(any(omit <- phi == 0)) {
		f <- y
		f[omit & y==mu] <- Inf
		f[omit & y!=mu] <- 0
		if(any(!omit))
			f[!omit] <- Recall(y[!omit], mu[!omit], phi[!omit], p)
		return(f)
	}
	if(p == 1)
		return(dpois(y/phi, lambda))
	if(p == 2)
		return(dgamma(y/beta, 1/phi)/beta)

	dinterp <- (2-p) * dpois( y/beta*(p-1), lambda*(2-p)) / beta * (p-1) +
	           (p-1) * dgamma(y/beta*(p-1), lambda*(2-p)) / beta * (p-1)

	omitn <- y < 0
	omit0 <- (y == 0 & !omitn)
	omitl <- (lambda > 100000 & !omitn & !omit0)
	omit1 <- (dinterp == 0 & !omitn & !omit0 & !omitl)
	if(any(omit <- omitn | omit0 | omitl | omit1)) {
		f <- y
		f[omitn] <- 0
		f[omit0] <- exp( - lambda[omit0])
		f[omitl] <- dinterp[omitl]
		f[omit1] <- 0
		if(any(!omit))
			f[!omit] <- Recall(y[!omit], mu[!omit], phi[!omit], p)
		return(f)
	}
#
#  Poisson mixture of gamma functions
#
	y <- y/beta
	k0 <- floor(min(lambda))
	if(k0 < 10)
		k0 <- 1
	k <- k0
	poi0 <- poi <- dpois(k, lambda)
	f <- gampoi <- gampoi0 <- dgamma(y, k * alpha)/beta * poi
	while(any(gampoi > gampoi0*1e-012)) {
		k <- k + 1
		poi <- (poi * lambda)/k
		gampoi <- dgamma(y, k * alpha)/beta * poi
		f <- f + gampoi
	}
	if(k0 == 1)
		return(f)
	k <- k0 - 1
	poi <- (poi0 * k0)/lambda
	gampoi <- dgamma(y, k * alpha)/beta * poi
	f <- f + gampoi
	while(k > 1 & any(gampoi > gampoi0*1e-012)) {
		poi <- (poi * k)/lambda
		k <- k - 1
		gampoi <- dgamma(y, k * alpha)/beta * poi
		f <- f + gampoi
	}
	f
}

ppoisgam <- function(y, mu = stop("No mean specified"), phi = 1, p = 1.5)
{
#  Cumulative distribution function for the Tweedie family with 1 <= p <= 2.
#  Gordon Smyth, 3 Mar 96, 13 Dec 97, 6 Dec 98
#
#  Check input arguments are admissable
#
	if(any(mu < 0)) stop("Mean of Tweedie family must be non-negative")
	if(any(phi < 0))
		stop("Dispersion of Tweedie family must be non-negative")
	if(length(p) > 1)
		stop("Index of Tweedie family must be scalar")
	if(p < 1 | p > 2) stop("Index of Tweedie family must be between 1 and 2")
#
#  Make sure arguments are the same length
#
	ny <- length(y)
	nmu <- length(mu)
	nphi <- length(phi)
	n <- max(c(ny, nmu, nphi))
	if(ny < n) y <- rep(y, length = n)
	if(nmu < n) mu <- rep(mu, length = n)
	if(nphi < n) phi <- rep(phi, length = n)
#
#  The Tweedie distribution can be expressed as a Poisson(lambda) sum of
#  gamma(alpha,beta)
#
	lambda <- mu^(2 - p)/(2 - p)/phi
	alpha <- (2 - p)/(p - 1)
	beta <- phi * (p - 1) * mu^(p - 1)
#
#  Special cases
#
	if(any(omit <- phi == 0)) {
		f <- y
		f[omit] <- y >= mu
		if(any(!omit))
			f[!omit] <- Recall(y[!omit], mu[!omit], phi[!omit], p)
		return(f)
	}
	if(p == 1)
		return(ppois(y/phi, lambda))
	if(p == 2)
		return(pgamma(y/beta, 1/phi))

	pinterp <- (2-p) * ppoiscc(y/beta*(p-1), lambda*(2-p)) +
	           (p-1) * pgamma( y/beta*(p-1), lambda*(2-p))

	omitn <- y < 0
	omit0 <- (y == 0 & !omitn)
	omitl <- (lambda > 100000 & !omitn & !omit0)
	omit1 <- (pinterp == 1 & !omitn & !omit0 & !omitl)
	if(any(omit <- omitn | omit0 | omitl | omit1)) {
		f <- y
		f[omitn] <- 0
		f[omit0] <- exp( - lambda[omit0])
		f[omitl] <- pinterp[omitl]
		f[omit1] <- 1
		if(any(!omit))
			f[!omit] <- Recall(y[!omit], mu[!omit], phi[!omit], p)
		return(f)
	}
#
#  Poisson mixture of gamma functions
#
	y <- y/beta
	k0 <- floor(min(lambda))
	if(k0 < 10)
		k0 <- 1
	k <- k0
	poi0 <- poi <- dpois(k, lambda)
	gampoi <- pgamma(y, k * alpha) * poi
	f <- exp( - lambda) + gampoi
	while(max(gampoi) > 1e-012) {
		k <- k + 1
		poi <- (poi * lambda)/k
		gampoi <- pgamma(y, k * alpha) * poi
		f <- f + gampoi
	}
	if(k0 == 1)
		return(f)
	k <- k0 - 1
	poi <- (poi0 * k0)/lambda
	gampoi <- pgamma(y, k * alpha) * poi
	f <- f + gampoi
	while(k > 1 & max(gampoi) > 1e-012) {
		poi <- (poi * k)/lambda
		k <- k - 1
		gampoi <- pgamma(y, k * alpha) * poi
		f <- f + gampoi
	}
	f
}

ppoiscc <- function(q, lambda = stop("no lambda arg"))
{
#	Poisson distribution function with continuity correction
#	Gordon Smyth, 6 Dec 98
#
	fq <- floor(q)
	d <- q - fq
	ppois(fq,lambda) + (d - 0.5) * dpois(fq + (d>0.5), lambda)
}

rpoisgam <- function(n, mu = stop("No mean specified"), phi = 1, p = 1.5)
{
#  Random variables for the Tweedie family with 1 <= p <= 2.
#  Gordon Smyth, 13 Dec 2001
#
#  Check input arguments are admissable
#
	if(length(n) > 1) n <- length(n)
	n <- floor(n)
	if(n < 1) return(NULL)
	if(any(mu < 0)) stop("Mean of Tweedie family must be non-negative")
	if(any(phi < 0))
		stop("Dispersion of Tweedie family must be non-negative")
	if(length(p) > 1) {
		p <- p[1]
		warning("p of length > 1. Only the first element used")
	}
	if(p < 1 | p > 2) stop("Index of Tweedie family must be between 1 and 2")
#
#  Make sure arguments are the same length
#
	nmu <- length(mu)
	nphi <- length(phi)
	if(nmu < n) mu <- rep(mu, length = n)
	if(nphi < n) phi <- rep(phi, length = n)
#
#  The Tweedie distribution can be expressed as a Poisson(lambda) sum of
#  gamma(alpha,beta)
#
	lambda <- mu^(2 - p)/(2 - p)/phi
	alpha <- (2 - p)/(p - 1)
	beta <- phi * (p - 1) * mu^(p - 1)
#
#  Special cases
#
	if(any(omit <- phi*mu^p == 0)) {
		r <- mu
		if(any(!omit))
			r[!omit] <- Recall(n-sum(omit), mu[!omit], phi[!omit], p)
		return(r)
	}
	if(p == 1)
		return(phi * rpois(n, lambda))
	if(p == 2)
		return(rgamma(n, shape=1/phi, rate=1/mu/phi))
#
#  Poisson mixture of gamma functions
#
	N <- rpois(n,lambda)
	r <- N
	if(any(N1 <- N > 0)) r[N1] <- rgamma(sum(N1), shape=N[N1]*alpha, rate=1/beta[N1])
	r
}
