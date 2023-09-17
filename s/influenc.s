influence <- function(object, ...)
UseMethod("influence")

influence.lm <- function(lm.obj)
{
#	Residuals and influence measures for lm objects
#	Useful also for glm and gam objects
#	GKS  Jan 96, 15 Oct 97, 27 Jan 98
#
	lms <- summary.lm(lm.obj)
	lmi <- lm.influence(lm.obj)
	h <- lmi$hat
	e <- residuals(lm.obj, type = "pearson")
	z <- e/sqrt(1 - h)/lms$sigma
	t <- e/sqrt(1 - h)/lmi$sigma
	k <- lms$df[1]
	d <- (z^2 * h)/k/(1 - h)
	b <- lms$coefficients[2:k, 1]
	se.b <- lms$coefficients[2:k, 2]
	bi <- lmi$coefficients[, 2:k]
	se.bi <- lmi$sigma %o% (se.b/lms$sigma)
	n <- sum(lms$df[1:2])
	dfbetas <- (rep(1, n) %o% b - bi)/se.bi
	data.frame(h, z, t, d, dfbetas)
}

plot.t <- function(inf)
{
#	Plot externally Studentized residuals from a linear model
#	GKS  Jan 1996
#
	plot(inf$t, type = "h", ylab = "Studentized Residuals")
	n <- length(inf$t)
	p <- dim(inf)[2] - 4
	t.crit <-  - qt(0.025/n, n - p - 2)
	lines(c(1, n), c(t.crit, t.crit), lty = 2)
	lines(c(1, n), c( - t.crit,  - t.crit), lty = 2)
}

plot.d <- function(inf)
#	Plot Cook's distances from a linear model
#	GKS  Jan 1996
#
	plot(inf$d, type = "h", ylab = "Cook's Distance")

