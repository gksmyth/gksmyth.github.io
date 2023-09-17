Digamma <- function(link = "log")
{
#	Digamma generalized linear model family
#	Gordon Smyth  3 July 1998
#
	name.link <- substitute(link)
	if(is.name(name.link))
		if(is.character(link))  # link is character variable
			name.link <- link
		else                    # link is name without quotes
			link <- name.link <- as.character(name.link)
	else
		if(is.call(name.link))  # power link
			name.link <- deparse(name.link)
	if(match(name.link, dimnames(glm.links)[[2]], F))
		link <- glm.links[, name.link]
	if(!is.null(link$name))
		name.link <- link$name

	var <- list(
		name = "2*[1/theta + trigamma(-theta)]",
		variance = varfun.digamma,
		deviance = function(mu, y, w = 1, residuals = F) {
			devi <- deviance.digamma(y,mu)
			if(residuals)
				sign(y - mu) * sqrt(abs(devi))
			else sum(devi)
		}
	)
	make.family("Digamma", link, var, name.link, "Trigamma")
}

kappa.digamma <- function(theta)
#	Cumulant function for the Digamma family
#	GKS  3 July 98
	2*( theta*(log(-theta)-1) + lgamma(-theta) )

meanval.digamma <- function(theta)
#	Mean value function for the Digamma family
#	GKS  3 July 98
	2*( log(-theta) - digamma(-theta) )

d2kappa.digamma <- function(theta)
#	2nd derivative of cumulant functio for Digamma family
#	GKS  3 July 98
	2*( 1/theta + trigamma(-theta) )

canonic.digamma <- function(mu) {
#	Canonical mapping for Digamma family
#	Solve meanval.digamma(theta) = mu for theta
#	GKS  3 July 98
#
#	Starting value from -log(-theta) =~ log(mu)
	mlmt <- log(mu)	
	theta <- -exp(-mlmt)

	for (i in 1:3) {
		mu1 <- meanval.digamma(theta)
		v <- d2kappa.digamma(theta)
		deriv <- -v/mu1*theta
		mlmt <- mlmt - log(mu1/mu)/deriv	
		theta <- -exp(-mlmt)
	}
	theta
}

varfun.digamma <- function(mu) {
#	Variance function for Digamma family
#	GKS  3 July 98
#
	theta <- canonic.digamma(mu)
	2*( 1/theta + trigamma(-theta) )
}

deviance.digamma <- function(y,mu) {
#	Unit deviance for Digamma family
#	GKS  3 July 98
#
	thetay <- canonic.digamma(y)
	theta <- canonic.digamma(mu)
	2*( y*(thetay-theta) - (kappa.digamma(thetay)-kappa.digamma(theta)) )
}

glm.weight <- function(link, variance)
{
#  This function fixes a bug in S-Plus 2000 Release 1.
#  It is not required in earlier or later versions of S-Plus.
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#  5 Nov 1999.
#
	default <- expression(w/((sqrt(family$variance(mu)) * family$deriv(mu))^2))
	dnames <- dimnames(glm.weights)
	if(!match(link, dnames[[1]], F))
		return(default)
	if(!match(variance, dnames[[2]], F))
		return(default)
	ww <- glm.weights[link, variance]
	if(as.character(ww) == "NULL")
		default
	else ww
}
