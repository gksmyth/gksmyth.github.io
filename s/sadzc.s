eppmsadzc <- function(lambda) {
# Saddle-point approximation based on the distribution with 6 nonzero cumulants
# Returns log-probability
# Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
# 24 May 1999.  Last revised 3 Nov 1999.

# Count lambda and test for special cases
n <- length(lambda)-1
if(n == 0) return(-lambda)
lambdan <- lambda[1:n]
if(any(lambdan <= 0)) return(-Inf)
if(any(is.na(lambda))) return(NA)
lambda1 <- lambda[n+1]
if(lambda1 < 0) lambda1 <- 0
l0 <- min(lambda)

# Find canonical parameter theta to make the tilted mean = 1
theta <- min(l0-1,mean(lambda)-n-1)
dif <- sum(1/(lambda-theta)) - 1
while (dif > 1e-15) {
	ddif <- sum(1/(lambda-theta)^2)
	theta <- theta - dif / ddif
	dif <- sum(1/(lambda-theta)) - 1
}

# Cumulants of tilted distribution
k <- rep(1,6)
k[2] <- sum(1/(lambda-theta)^2)
k[3] <- 2*sum(1/(lambda-theta)^3)
k[4] <- 6*sum(1/(lambda-theta)^4)
k[5] <- 24*sum(1/(lambda-theta)^5)
k[6] <- 120*sum(1/(lambda-theta)^6)

# Compute probability at 1
log(dzeroc(1,k))-theta+sum(log(lambdan/(lambdan-theta)))-log(lambda[n+1]-theta)
}
