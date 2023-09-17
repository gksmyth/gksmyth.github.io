digamma <- function(z)
{
#  Digamma function: first derivative of log(gamma(x)).
#  Original function by Bill Venables.
#  Modified by GKS to accepts NAs, 19 Jan 1998.
#
	if(any(omit <- is.na(z) | Re(z) <= 0)) {
		ps <- z
		ps[omit] <- NA
		if(any(!omit))
			ps[!omit] <- Recall(z[!omit])
		return(ps)
	}
	if(any(small <- Mod(z) < 5)) {
		ps <- z
		x <- z[small]
		ps[small] <- Recall(x + 5) - 1/x - 1/(x + 1) - 1/(x + 2) - 1/(x + 3) - 1/(x + 4)
		if(any(!small))
			ps[!small] <- Recall(z[!small])
		return(ps)
	}
	x <- 1/z^2
	tail <- ((x * (-1/12 + ((x * (1/120 + ((x * (-1/252 + ((
		x * (1/240 + ((x * (-1/132 + ((x * (691/32760 + (
		(x * (-1/12 + (3617 * x)/8160)))))))))))))))))))))
	log(z) - 1/(2 * z) + tail
}

trigamma <- function(z)
{
#  Trigamma function: second derivative of log(gamma(x)).
#  Original function by Bill Venables.
#  Modified by GKS to accepts NAs, 19 Jan 1998.
#
	if(any(omit <- is.na(z) | Re(z) <= 0)) {
		ps <- z
		ps[omit] <- NA
		if(any(!omit))
		ps[!omit] <- Recall(z[!omit])
		return(ps)
	}
	if(any(small <- Mod(z) < 5)) {
		ps <- z
		x <- z[small]
		ps[small] <- Recall(x + 5) + 1/x^2 + 1/(x + 1)^2 +
			1/(x + 2)^2 + 1/(x + 3)^2 + 1/(x + 4)^2
		if(any(!small))
		ps[!small] <- Recall(z[!small])
		return(ps)
	}
	x <- 1/z^2
	tail <- 1 + (x * (1/6 + (x * (-1/30 + (x * (1/42 + (x * (
		-1/30 + (x * (5/66 + (x * (-691/2370 + (x * (7/6 -
		(3617 * x)/510))))))))))))))
	1/(2 * z^2) + tail/z
}

logmdig <- function(z)
{
#  log(z) - digamma(z)
#  Saves computation of log(z) and avoids subtractive cancellation when z is large
#  GKS  19 Jan 98
#
	if(any(omit <- is.na(z) | Re(z) <= 0)) {
		ps <- z
		ps[omit] <- NA
		if(any(!omit))
			ps[!omit] <- Recall(z[!omit])
		return(ps)
	}
	if(any(small <- Mod(z) < 5)) {
		ps <- z
		x <- z[small]
		ps[small] <- log(x/(x+5)) + Recall(x + 5) + 1/x + 1/(x + 1) + 1/(x + 2) + 1/(x + 3) + 1/(x + 4)
		if(any(!small))
			ps[!small] <- Recall(z[!small])
		return(ps)
	}
	x <- 1/z^2
	tail <- ((x * (-1/12 + ((x * (1/120 + ((x * (-1/252 + ((
		x * (1/240 + ((x * (-1/132 + ((x * (691/32760 + (
		(x * (-1/12 + (3617 * x)/8160)))))))))))))))))))))
	1/(2 * z) - tail
}
