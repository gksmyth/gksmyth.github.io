dtweedie <- function(y, mu, phi=1, power=1.5)
{
#  Density or probability mass function for Tweedie
#  exponential dispersion models.
#
#  mu is the mean, and phi is the dispersion of the distribution
#  We have var(y) = phi * mu^p
#
#  GKS  3 June 98.  Last revised 10 Nov 99.
#
	if(power == 0) {
		return( dnorm(y, mu, phi) )
	}
	if(power == 1) {
		return( dpois(y/phi, mu/phi) )
	}
	if(power>1 & power<2) {
		return( dpoisgam(y, mu, phi, power) )
	}
	if(power == 2) {
		alpha <- 1/phi
		beta <- mu/alpha
		return( dgamma(y/beta, alpha)/beta )
	}
	if(power == 3) {
		return( dinvgauss(y, mu, 1/phi) )
	}
	saddle.tweedie(y, mu, phi, power)
}

ptweedie <- function(q, mu, phi=1, power=1.5)
{
#  Cumulative distribution function for Tweedie
#  exponential dispersion models.
#
#  mu is the mean, and phi is the dispersion of the distribution
#  We have var(y) = phi * mu^p
#
#  GKS  30 Mar 99
#
	if(power == 0) {
		return( pnorm(q, mu, phi) )
	}
	if(power == 1) {
		return( ppois(q/phi, mu/phi) )
	}
	if(power>1 & power<2) {
		return( ppoisgam(q, mu, phi, power) )
	}
	if(power == 2) {
		alpha <- 1/phi
		beta <- mu/alpha
		return( pgamma(q/beta, alpha) )
	}
	if(power == 3) {
		return( pinvgauss(q, mu, 1/phi) )
	}
#
#	Modified residual r* tail approximation
	tail.tweedie(q, mu, phi, power)
}

qtweedie <- function(p, mu, phi=1, power=1.5)
{
#  Quantiles for Tweedie exponential dispersion models
#
#  mu is the mean and phi is the dispersion of the distribution
#  We have var(y) = phi * mu^power
#
#  This function is currently reasonable only for the normal, Poisson, gamma
#  and inverse-Gaussian special cases.
#
#  GKS  28 Mar 99, 4 Oct 99
#
	if(power == 0) {
		return( qnorm(p, mu, phi) )
	}
	if(power == 1) {
		return( phi * qpois(p, mu/phi) )
	}
	if(power == 2) {
		alpha <- 1/phi
		beta <- mu/alpha
		return( beta * qgamma(p, alpha) )
	}
	if(power == 3) {
		return( qinvgauss(p, mu, 1/phi) )
	}
#
#	Simple normal approximation - need to improve!
	qnorm(p, mu, phi*mu^power)
}

saddle.tweedie <- function(y, mu, phi, p)
{
#  Saddlepoint approximation for Tweedie densities
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#  3 June 1998.  Last revised 10 Nov 1999.
#
	y1 <- ifelse(y>0,y,1/6)
	exp( - devi.tweedie(y, mu, p)/2/phi)/sqrt(2 * pi * phi * y1^p)
}

devi.tweedie <- function(y, mu = 1, p = 0)
{
#  Unit deviance of Tweedie family
#  GKS  3 June 98
#
	y1 <- y + (y == 0)
	if(p == 1)
		theta <- log(y1/mu)
	else theta <- (y1^(1 - p) - mu^(1 - p))/(1 - p)
	if(p == 2)
		kappa <- log(y1/mu)
	else kappa <- (y^(2 - p) - mu^(2 - p))/(2 - p)
	2 * (y * theta - kappa)
}

tail.tweedie <- function(y, mu = 1, phi = 1, p = 0)
{
#  Tail area approximation based on modified deviance residual
#  See Jørgensen (1997) page 143
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#  28 April 2000
#
	if(phi==0) return(NA)
	e <- y-mu
	s <- sqrt(phi)
	z <- y/mu
	u <- mu^(1-p/2) * z^(p/2) * ifelse(p==1,log(z),(z^(1-p)-1)/(1-p))
	r <- sign(e) * sqrt(devi.tweedie(y,mu,p))
	rmod <- r/s + ifelse(abs(e/mu) > 1e-8, s/r*log(u/r), p*s*mu^(p/2-1)/6)
	pnorm(rmod)
}
