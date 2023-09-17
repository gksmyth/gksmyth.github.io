sym.plot <- function(y)
{
# Symmetry plot of a numeric vector
# Gordon Smyth  Jan 1996
#
   y <- sort(y)
   m <- median(y)
   p <- (length(y)+1)/2
   bot <- m - y[ floor(p-0.5):1 ]
   top <- y[ ceiling(p+0.5):length(y) ] - m
   plot(bot,top,xlab="Dist below median",ylab="Dist above median")
   a <- max(min(bot),min(top))
   b <- min(max(bot),max(top))
   lines(c(a,b),c(a,b))
}
