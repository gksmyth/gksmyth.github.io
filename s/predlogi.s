pointwise.logit <- function(glm.obj,newdata,coverage=0.99)
{
# Confidence intervals for prediction for logistic regression
# Gordon Smyth  Feb 1996
#
	lp <- predict(glm.obj,newdata,se=T)
	a <- (1+coverage)/2
	margin.error <- lp$se.fit / lp$residual.scale * qnorm(a)
	lp.upper <- lp$fit + margin.error
	lp.lower <- lp$fit - margin.error
	upper <- exp(lp.upper)
	upper <- upper/(1+upper)
	lower <- exp(lp.lower)
	lower <- lower/(1+lower)
	fit <- exp(lp$fit)
	fit <- fit/(1+fit)
	list(upper=upper,fit=fit,lower=lower)
}
