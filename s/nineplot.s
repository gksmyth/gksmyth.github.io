nineplot <- function(y, ...)
{
# NINEPLOT
# Normal probability plot surrounded by reference plots
# GKS  Jul 96, 29 Dec 96
#
	oldpar <- par(mfrow = c(3, 3), mar = rep(1, 4), col = 1)
	on.exit(par(oldpar))
	doqqplot <- function(x, p, ...)
	{
		q <- qnorm(p)[order(order(x))]
		.S(plot(range(q), range(x), type = "n", axes = F, ...), "plot")
		box(...)
		points(q, x, ...)
		qqline(x, ...)
	}
	y <- y[!is.na(y)]
	n <- length(y)
	p <- ppoints(n)
	for(i in 1:4) {
		z <- rnorm(n)
		doqqplot(z, p, col = 1)
	}
	doqqplot(y, p, col = 2)
	for(i in 1:4) {
		z <- rnorm(n)
		doqqplot(z, p, col = 1)
	}
	invisible()
}

