# Weibull Hazard Plot
# Developed by Dr. Reuel Smith, 2023

plot.hazard.wbl <- function(xi, rc = NULL, pp = "Blom", confid = 0.95, sided = "twosided", xlabel1 = "X"){
  library(pracma)
  library(mvtnorm)
  library(ggplot2)

  # Get nonparametric estimates
  x_R <- plotposit.select(xi, rc, pp)

  # Run Weibull fit
  params <- fit.wbl(xi, rc, pp, confid,sided)
  # Run histogram analysis to get plotting range
  histoutput <- hist(rweibull(10000, params$wblparams[2], params$wblparams[1]),1000)
  # Pull the x and y values
  xrange <- linspace(min(histoutput$breaks),max(histoutput$breaks),1000)
  yhx <- (params$wblparams[2]/params$wblparams[1])*((xrange/params$wblparams[1])^(params$wblparams[2]-1))
  # yhx <- dweibull(xrange, params$wblparams[2], params$wblparams[1])/(1 - pweibull(xrange, params$wblparams[2], params$wblparams[1]))
  # Confidence bounds calculation by transformation
  u_transform <- params$wblparams[2]*(log(xrange) - log(params$wblparams[1]))
  var_u <- ((params$wblparams[2]/params$wblparams[1])^2)*params$varcovmat[1,1] +
    ((u_transform/params$wblparams[2])^2)*params$varcovmat[2,2] +
    2*(u_transform/params$wblparams[1])*params$varcovmat[1,2]
  crit <- qnorm((1 + confid)/2)
  u_low<-u_transform - (crit * sqrt(var_u))
  u_high<-u_transform + (crit * sqrt(var_u))
  yhx_low <- exp(-exp(u_low))
  yhx_high <- exp(-exp(u_high))

  # Form data frame
  df <- data.frame(X = xrange, YHaz = yhx, YHazlow = yhx_low, YHazhigh = yhx_high, best_fit = rep("Fitted",1000))
  # return(list(df))
  df2 <- data.frame(Xpts = x_R[,1][1:(length(x_R[,1])-1)], Ypts = x_R[,4][1:(length(x_R[,1])-1)], nonparametric_data = rep("Failure Data",length(x_R[,1])-1))
  plotout<-ggplot() +
    geom_line(data=df, aes(X,YHaz), colour = 'blue', size = 0.9, linetype = "dashed") +
    xlab(xlabel1) +
    ylab("Hazard Rate")
  # plotout <- plotout + geom_ribbon(data=df, aes(ymin=YHazlow, ymax=YHazhigh, x=X), alpha=0.5, fill = "blue")
  plotout <- plotout + geom_point(data=df2, aes(Xpts,Ypts, shape = nonparametric_data), colour = 'red', size = 2.2)
  return(list(wbloutput = params,Hazardplot = plotout))
}
