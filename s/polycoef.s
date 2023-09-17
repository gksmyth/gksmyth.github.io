polycoef <- function(x) {
#   Produces a vector whose elements are the coefficients of
#	the polynomial whose roots are the elements of x.
#	polyroot and polycoef are inverse functions of each other,
#	up to ordering, scaling, and roundoff error.

#	Gordon Smyth, gks@maths.uq.edu.au, 4 June 1999

# Strip out infinities
e <- x[ is.finite(x) ]

# Expand recursion formula
n <- length(e)
a <- c(rep(0,n),1)
for (j in n:1) a[j:n] <- a[j:n] - e[j]*a[(j+1):(n+1)]

# The result should be real if the roots are complex conjugates
pos <- sort(e[Im(e)>0])
neg <- sort(Conj(e[Im(e)<0]))
if(length(pos)==length(neg))
	if(all(pos==neg))
		a <- Re(a)

a
}
