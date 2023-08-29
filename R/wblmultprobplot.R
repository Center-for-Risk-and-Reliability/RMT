# Multiple Weibull Probability Plot
# Developed by Dr. Reuel Smith, 2021-2022

multiprobplot.wbl <- function(data,pp) {
  databystress<-checkstress(data)

  sort.xircdata(data)
  # Data pull
  XB<-xiRFblock[,1]
  FB <- log(log(1/xiRFblock[,3]))

  # Initialize Percent Failure ticks
  Pticks1 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,c(1:10),10*c(2:9),95,99,99.9)
  Pticks <- log(log(1/(1-Pticks1/100)))
  Pticks1label <- c(0.1,0.2,0.3,"",0.5,"","","","",1,2,3,"",5,"","","","",10*c(1:9),95,99,99.9)

  # For now just do a single stress plot
  # Note for help file: xiFR is to come from a plotposit series function
  ttfc<-probplotparam.wbl(xiRFblock[,1],xiRFblock[,3])
  # I will add a multi-stress example a little later

  # Computes the upper and lower bound for the TTF axis in terms of log-time
  signs1 <- c(floor(log10(min(ttfc[[1]]))):ceiling(log10(max(ttfc[[1]]))))
  logtimes1 <- 10^signs1
  Pticks1X <- c(1:(9*length(logtimes1)-8))
  Pticks1X[1] <- logtimes1[1]
  Pticks1Xlabel <- Pticks1X

  for(i2 in 1:(length(signs1)-1)){
    Pticks1X[(9*i2-7):(9*(i2+1)-8)] <- logtimes1[i2]*c(2:10)
    Pticks1Xlabel[(9*i2-7):(9*(i2+1)-8)] <- c("","","","","","","","",logtimes1[i2+1])
  }

  # Plot
  plot(XB, FB, log="x", col="blue",
       xlab="Time (hours)", ylab="Percent Failure", pch=16,
       xlim=c(10^min(signs1), 10^max(signs1)), ylim=c(min(ttfc[[2]]), max(ttfc[[2]])), axes=FALSE)
  lines(ttfc[[1]], ttfc[[2]], col="blue")
  axis(1, at=Pticks1X,labels=Pticks1Xlabel, las=2, cex.axis=0.7)
  axis(2, at=Pticks,labels=Pticks1label, las=2, cex.axis=0.7)

  #Add horizontal and vertical grid
  abline(h = Pticks, lty = 2, col = "grey")
  abline(v = Pticks1X,  lty = 2, col = "grey")

  legend(10^signs1[1], log(log(1/(1-.99))), legend=c("Data", "Weibull best fit line"),
         col=c("blue", "blue"), pch=c(16,-1), lty=c(0,1), cex=0.8,
         text.font=4, bg='white')
  return(ttfc[[3]])
}
