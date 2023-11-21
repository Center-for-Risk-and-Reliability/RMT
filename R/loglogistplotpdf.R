# Log-logistic PDF Plot
# Developed by Dr. Reuel Smith, 2023

plot.pdf.loglogist <- function(xi, rc = NULL, pp = "Blom", confid = 0.95, sided = "twosided", xlabel1 = "X"){
  library(pracma)
  library(ggplot2)

  # Run Log-logistic fit
  params <- fit.loglogist(xi, rc, pp, confid,sided)
  # Run histogram analysis to get plotting range
  histoutput <- hist(exp(params$loglogistparams[1] + params$loglogistparams[2]*log(((1/runif(10000,0,1)) - 1)^-1)),1000)
  # Pull the x and y values
  xrange <- linspace(min(histoutput$breaks),max(histoutput$breaks),1000)
  ypdf <- (1/xrange)*dlogis(log(xrange), params$loglogistparams[1], params$loglogistparams[2])

  # Form data frame
  df <- data.frame(X = xrange, YPDF = ypdf, data = rep("Fitted",1000))
  plotout<-ggplot() +
    geom_line(data=df, aes(X,YPDF, colour = data), colour = 'blue', size = 0.9, linetype = "dashed") +
    xlab(xlabel1) +
    ylab("Probability Density")
  return(list(loglogistoutput = params,pdfplot = plotout))
}
