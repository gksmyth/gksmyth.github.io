eppmsadga <- function(lambda,second=T) {
# Saddle-point approximation based on gamma distribution
# Returns log-probability
# Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
# 3 November 1999

# Count lambda and test for Poisson case
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
iter <- 0
while (dif > 1e-13 && iter < 10) {
	iter <- iter + 1
	ddif <- sum(1/(lambda-theta)^2)
	theta <- theta - dif / ddif
	dif <- sum(1/(lambda-theta)) - 1
}

# Compute probability at 1
k1 <- sum(1/(lambda-theta))
k2 <- sum(1/(lambda-theta)^2)
b <- k2
a <- k1/b
ans <- (a-1)*log(k1)-k1/b-a*log(b)-lgamma(a)-k1*theta+sum(log(lambdan/(lambdan-theta)))-log(lambda[n+1]-theta)

# Second term correction
if(second) {
	k3 <- sum(1/(lambda-theta)^3)
	k4 <- sum(1/(lambda-theta)^4)
	K3 <- k2*b
	K4 <- K3*b
	ans <- ans + 1/8*(k4-K4)/k2^2 - 5/24*(k3^2-K3^2)/k2^3
}
ans
}
