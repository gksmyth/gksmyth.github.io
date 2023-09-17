tweedie <- function(var.power = 0, link.power = 1 - var.power)
{
#	Tweedie generalized linear model family
#	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	March 1996.  Revised 6 Aug 2001.
#
	if (link.power != 0) tw.lnk <- list(
		names = "Box-Cox: mu^q",
		link = substitute(
			function(mu, q = l.p) mu^q,
			list(l.p = link.power)),
		inverse = substitute(
			function(eta, q = l.p) eta^(1/q),
			list(l.p = link.power)),
		deriv = substitute(
			function(mu, q = l.p) q*mu^(q-1),
			list(l.p = link.power)),
		initialize = expression({
			y <- as.numeric(y)
			mu <- y + (y == 0)/6}))

	if (link.power == 0) tw.lnk <- list(
		names = "Log: log(mu)",
		link = function(mu) log(mu),
		inverse = function(eta) care.exp(eta),
		deriv = function(mu) 1/mu,
		initialize = expression({
			y <- as.numeric(y)
			mu <- y + (y == 0)/6}))

	tw.var <- list(
		names = "mu^p",
		variance = substitute(
			function(mu, p = v.p) mu^p,
			list(v.p = var.power)),
		deviance = substitute(
			function(mu, y, w, residuals = F, p = v.p) {
				y1 <- y + (y == 0)
				if (p == 1)
					theta <- log(y1/mu)
				else
					theta <- ( y1^(1-p) - mu^(1-p) ) / (1-p)
				if (p == 2)
					kappa <- log(y1/mu)
				else
					kappa <- ( y^(2-p) - mu^(2-p) ) / (2-p)
				devi <- y*theta - kappa
				if(residuals) sign(y-mu)*sqrt(2*w*abs(devi)) else 2*sum(w*devi)},
			list(v.p = var.power)))

	make.family("Tweedie", link = tw.lnk, variance = tw.var)
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
