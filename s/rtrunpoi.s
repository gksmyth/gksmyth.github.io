rtrunpois <- function(n, lambda = stop("no lambda arg"), origin=0) {
#  Random deviates from truncated Poisson distribution
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#  28 May 2000

if(length(n)>1) n <- length(n)
n <- floor(n)
origin <- floor(origin)

if(origin <= 0) return(rpois(n,lambda))
if(n <= 0) return(NULL)

p <- ppois(origin-1,lambda)
y <- rpois(ceiling(n/(1-p)),lambda)
y <- y[y>=origin]

n1 <- length(y)
if(n1 >= n)
	return( y[1:n] )
else
	return( c(y, Recall(n-n1,lambda,origin) ))
}
