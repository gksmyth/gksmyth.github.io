plotcircles <- function(x, y, size = 1, xlab = NULL, ylab = NULL, ...)
{
#	Scatterplot with circle size indexing third variable
#	GKS 10 September 1997
#
	call <- match.call()
	if(is.null(xlab))
		xlab <- deparse(call$x)
	if(is.null(ylab))
		ylab <- deparse(call$y)
	plot(x, y, xlab = xlab, ylab = ylab, type = "n", ...)
	symbols(x, y, add = T, circles = sqrt(size), inches = 0.25, ...)
	invisible()
}

