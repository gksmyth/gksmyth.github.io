eppmsadnb <- function(lambda,second=T) {
#  Saddle-point approximation based on negative binomial distribution
#  Returns log-probability
#  Gordon Smyth, U of Queensland, gks@maths.uq.edu.au
#  24 May 1999.  Last revised 3 Nov 1999.

# Count lambda and test for Poisson case
n <- length(lambda)-1
if(n == 0) return(-lambda)
lambdan <- lambda[1:n]
if(any(lambdan <= 0)) return(-Inf)
if(any(is.na(lambda))) return(NA)
lambda1 <- lambda[n+1]
if(lambda1 < 0) lambda1 <- 0
l0 <- min(lambda)
ln <- max(lambda)
if(ln-l0 < 1e-8) return(log(dpois(n,mean(lambda))))

# Find canonical parameter theta to make the tilted mean = 1
theta <- min(l0-1,mean(lambda)-n-1)
dif <- sum(1/(lambda-theta)) - 1
while (dif > 1e-15) {
	ddif <- sum(1/(lambda-theta)^2)
	theta <- theta - dif / ddif
	dif <- sum(1/(lambda-theta)) - 1
}

# Find negative binomial parameters to match waiting time mean and variance
k1 <- sum(1/(lambda-theta))
k2 <- sum(1/(lambda-theta)^2)
a <- (ln-l0)/n
b <- (l0-theta)/a
sum1 <- sum(1/(b+0:n))
sum2 <- sum(1/(b+0:n)^2)
dif <- k1^2*sum2/sum1^2-k2
tol <- (1e-14)*k2
while (dif < -tol) {
	b <- b/2
	sum1 <- sum(1/(b+0:n))
	sum2 <- sum(1/(b+0:n)^2)
	dif <- k1^2*sum2/sum1^2-k2
}
while (dif > tol) {
	sum3 <- sum(1/(b+0:n)^3)
	ddif <- 2*k1^2/sum1^2*(sum2^2/sum1-sum3)
	b <- b - dif / ddif
	sum1 <- sum(1/(b+0:n))
	sum2 <- sum(1/(b+0:n)^2)
	dif <- k1^2*sum2/sum1^2-k2
}

# Compute probability at 1
a <- sum1
if(a < 1e-5)
	expa <- k1*(1-a*k1/2+(a*k1)^2/6)
else
	expa <- (1-exp(-a*k1))/a
ans <- -a*b*k1-lgamma(n+1)+sum(log(a*(b+0:n)))+n*log(expa)-
	theta+sum(log(lambdan/(lambdan-theta)))-log(lambda[n+1]-theta)

# Second term correction
if(second) {
	k3 <- sum(1/(lambda-theta)^3)
	k4 <- sum(1/(lambda-theta)^4)
	K3 <- sum(1/(a*(b+0:n))^3)
	K4 <- sum(1/(a*(b+0:n))^4)
	ans <- ans + 3/4*(k4-K4)/k2^2 - 5/6*(k3^2-K3^2)/k2^3
}
ans
}
