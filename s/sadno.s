eppmsadno <- function(lambda,second=T) {
# Saddle-point approximation based on normal distribution
# Returns log-probability
# Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
# 19 May 1999.  Last revised 3 Nov 1999.

# Count lambda and test for Poisson case
n <- length(lambda)-1
if(n > 0) {
	lambdan <- lambda[1:n]
	if(any(lambdan <= 0)) return(-Inf)
}
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

# Compute probability at 1
k2 <- sum(1/(lambda-theta)^2)
ans <- -0.5*log(2*pi*k2)-theta-log(lambda[n+1]-theta)
if(n>0) ans <- ans+sum(log(lambdan/(lambdan-theta)))

# Second term correction
if(second) {
	k3 <- sum(1/(lambda-theta)^3)
	k4 <- sum(1/(lambda-theta)^4)
	rho3 <- k3/k2^1.5
	rho4 <- k4/k2^2
	ans <- ans + 3/4*rho4 - 5/6*rho3^2
}
ans
}
