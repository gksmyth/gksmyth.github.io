add.var <- function(lm.obj,...)
{
# Added variable plots for multiple regression
# Gordon Smyth  Feb 1996
#
	b <- coefficients(lm.obj)
	y <- formula(lm.obj)[[2]]
	k <- drop1(lm.obj,keep=c("x.residuals","residuals"),...)
	for (x in dimnames(k$keep)[[1]]) {
		plot(k$keep[[x,"x.resid"]],k$keep[[x,"resid"]],
			xlab=paste("Adjusted",x),
			ylab=paste("Adjusted",y),
			main=paste("Added Variable Plot for",x))
			abline(0,b[x])}
}
