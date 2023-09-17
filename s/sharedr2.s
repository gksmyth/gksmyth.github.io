shared.r2 <- function(x)
{
#	Shared R-squared (a measure of collinearity)
#	Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#	Feb 1996.
#
	n <- dim(x)[1]
	x <- x - rep(1,n) %o% apply(x,2,mean)
	xx <- t(x) %*% x
	xxinv <- solve(xx)
	tss <- diag(xx)
	varb <- diag(xxinv)
	r2 <- 1 - 1/(tss*varb)
	names(r2) <- dimnames(x)[[2]]
	r2
}
