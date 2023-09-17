polyval <- function(x,a) {
#   Evaluates at  X  the polynomial with coefficients given
#	by the elements of  A .  If  N = length(A)  then
#	POLYVAL <- A[1] + A[2]*X + ... + A[N]*X^(N-1)
#
#	Gordon Smyth, gks@maths.uq.edu.au, 6 June 1999
#
n <- length(a)
if(n==0) return(NULL)
y <- a[n]
if(n==1) return(y)
for (j in (n-1):1) y <- a[j] + x*y
y
}
