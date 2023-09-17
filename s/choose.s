choose <- function(n, x)
{
	exp(lgamma(n + 1) - lgamma(x) - lgamma(n - x))
}
