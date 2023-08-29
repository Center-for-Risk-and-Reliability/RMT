# Weibull PDF Plot
# Developed by Dr. Reuel Smith, 2023

plot.pdf.wbl <- function(xi, rc = NULL, pp = "Blom", confid = 0.95, sided = "twosided", xlabel1 = "X"){
  library(pracma)
  library(ggplot2)

  # Run Weibull fit
  params <- fit.wbl(xi, rc, pp, confid,sided)
  # Run histogram analysis to get plotting range
  histoutput <- hist(rweibull(10000, params$wblparams[2], params$wblparams[1]),1000)
  # Pull the x and y values
  xrange <- linspace(min(histoutput$breaks),max(histoutput$breaks),1000)
  ypdf <- dweibull(xrange, params$wblparams[2], params$wblparams[1])

  # Form data frame
  df <- data.frame(X = xrange, YPDF = ypdf, data = rep("Fitted",1000))
  plotout<-ggplot() +
    geom_line(data=df, aes(X,YPDF, colour = data), colour = 'blue', size = 0.9, linetype = "dashed") +
    xlab(xlabel1) +
    ylab("Probability Density")
  return(list(wbloutput = params,pdfplot = plotout))
}
