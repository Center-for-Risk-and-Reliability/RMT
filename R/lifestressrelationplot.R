# Life-Stress Relationship Plot Generator
# Developed by Dr. Reuel Smith, 2021-2022

lifestress.relationplot <- function(data,ls,dist,pp,Smin,Smax,SUse,therm,confid,Llab,Slab) {
  #Load plotly library for 3D plotting
  library(plotly)

  # Check first that the data has multiple accelerated stress levels
  if(length(pp)<3) {
    stop('Need more than one stress level to generate relationship plot.')
  }

  # Check to see that the vectors Smin, Smax, and SUse are the same length to determine
  # if life is based on one stress or multiple stress types
  if(isTRUE(length(Smin)==length(Smax)&&length(Smin)==length(SUse))){
    S_no<-length(Smax)
  } else {
    stop('Make sure that your Smin, Smax and SUse inputs are the same length.')
  }

  Sminmaxcheck<-rep(0,S_no)

  for (i2 in 1:S_no){
    if(isTRUE(Smax[i2]>Smin[i2])){
      Sminmaxcheck[i2]<-0
    } else{
      Sminmaxcheck[i2]<-1
    }
  }

  # Check to see if the order of stress max and min is correct (this for one stress)
  if(sum(Sminmaxcheck)>0) {
    stop('Your S_max is less than your S_min.')
  }

  # Check to see if a confidence bound was entered.  95% is the default.
  if(missing(confid)){
    confid <- 0.95
  }

  # Sort the stress data for MLE step
  xircSxiSrc <- sort.xircstressdata(data)

  # Compute the LSQ and MLE data
  LSQoutput <- lifestress.LSQest(ls,dist,pp)
  MLEoutput <- lifestress.MLEest(LSQoutput[[3]],ls,dist,xircSxiSrc[[1]],xircSxiSrc[[3]],xircSxiSrc[[2]],xircSxiSrc[[4]],confid)

  # ==========================================================================
  # Here is where we generate the relationship plots
  # ==========================================================================
  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential") {
    ishift<-0
  } else {
    ishift<-1
  }

  # Setup or select the Life-Stress function based on parameter estimates for theta and your
  # minimum and maximum stress values
  if (length(Smax)==1){
    Sline<-seq(Smin,Smax,(Smax-Smin)/100)
  }
  if (length(Smax)==2){
    S1set<-seq(Smin[1],Smax[1],(Smax[1]-Smin[1])/100)
    S2set<-seq(Smin[2],Smax[2],(Smax[2]-Smin[2])/100)
    S1line <- S1set%o%rep(1,101)
    S2line <- rep(1,101)%o%S2set
    #S1line <- rep(S1set, each = 100, len = 10100)
    #S2line <- rep(S2set,100)
    #Sline<-matrix(c(S1line,S2line),nrow=10100,ncol=2,byrow = FALSE)
  }

  #
  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    life <- function(theta,SF) {
      theta[ishift+2] + SF*theta[ishift+1]
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Linear Relationship Plot Stress vs. Life
    plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=c("Characteristic Stress (",Slab,")"),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
    lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(SUse,UselifeLSQ, pch=17, col="blue")
    points(SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b

    life <- function(theta,SF) {
      theta[ishift+2]*exp(SF*theta[ishift+1])
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Exponential Relationship Plot Stress vs. Life
    plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=c("Characteristic Stress (",Slab,")"),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
    lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(SUse,UselifeLSQ, pch=17, col="blue")
    points(SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HaS to be in Kelvin for this to work
    K<-8.617385e-5
    life <- function(theta,SF) {
      theta[ishift+2]*exp(theta[ishift+1]/(K*SF))
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Arrhenius Relationship Plot Inverse Stress vs. Life
    plot(1/Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=expression(paste("Characteristic Stress"^"-1","(","K"^"-1",")")),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(1/Smax,1/Smin),axes=FALSE)
    lines(1/Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(1/LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(1/SUse,UselifeLSQ, pch=17, col="blue")
    points(1/SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topleft", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    life <- function(theta,SF) {
      (theta[ishift+2]/SF)*exp(theta[ishift+1]/SF)
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Eyring Relationship Plot Inverse Stress vs. Life
    plot(1/Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=expression(paste("Characteristic Stress"^"-1","(","K"^"-1",")")),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(1/Smax,1/Smin),axes=FALSE)
    lines(1/Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(1/LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(1/SUse,UselifeLSQ, pch=17, col="blue")
    points(1/SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topleft", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    life <- function(theta,SF) {
      (1/SF)*exp(-(theta[ishift+1] - (theta[ishift+2]/SF)))
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Eyring 2 Relationship Plot Inverse Stress vs. Life
    plot(1/Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=expression(paste("Characteristic Stress"^"-1","(","K"^"-1",")")),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(1/Smax,1/Smin),axes=FALSE)
    lines(1/Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(1/LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(1/SUse,UselifeLSQ, pch=17, col="blue")
    points(1/SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topleft", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    life <- function(theta,SF) {
      theta[ishift+2]*(SF^theta[ishift+1])
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Power Relationship Plot Stress vs. Life
    plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=c("Characteristic Stress (",Slab,")"),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
    lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(SUse,UselifeLSQ, pch=17, col="blue")
    points(SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    life <- function(theta,SF) {
      theta[ishift+2]*(SF^-theta[ishift+1])
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Inverse Power Relationship Plot Stress vs. Life
    plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=c("Characteristic Stress (",Slab,")"),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
    lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(SUse,UselifeLSQ, pch=17, col="blue")
    points(SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    life <- function(theta,SF) {
      theta[ishift+1]*log(SF) + theta[ishift+2]
    }
    UselifeLSQ <- life(LSQoutput[[3]],SUse)
    UselifeMLE <- life(MLEoutput[[1]],SUse)

    # Logarithmic Relationship Plot Stress vs. Life
    plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
         xlab=c("Characteristic Stress (",Slab,")"),
         ylab=c("Characteristic Life (",Llab,")"),
         type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
    lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
    points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
    points(SUse,UselifeLSQ, pch=17, col="blue")
    points(SUse,UselifeMLE, pch=18, col="black")
    axis(1, at=NULL, las=2, cex.axis=0.7)
    axis(2, at=NULL, las=2, cex.axis=0.7)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
    legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
           col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
           text.font=4, bg='white')
  }

  if (ls=="MultiStress") {
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    life <- function(theta,SF) {
      exp(theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF))
    }
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    lifeUse <- function(theta,SF) {
      theta[ishift+1]*exp((theta[ishift+2]/SF[1]) + (theta[ishift+3]/SF[2]))
    } # Single stress entry life (for Use Life primarily)
    life <- function(theta,SF1,SF2) {
      theta[ishift+1]*exp((theta[ishift+2]/SF1) + (theta[ishift+3]/SF2))
    } # Multiple stress entry life (for vectors)
    UselifeLSQ <- lifeUse(LSQoutput[[3]],SUse)
    UselifeMLE <- lifeUse(MLEoutput[[1]],SUse)

    # Temperature-Humidity Relationship Plot Stress vs. Life
    persp(x=S1set,y=S2set,z=life(MLEoutput[[1]],S1line,S2line),theta = 100, phi = 50,col="cyan",
          shade=0.4, nticks=5, ticktype="detailed",
          xlab=c("Characteristic Temperature (",Slab[1],")"),
          ylab=c("Relative Humidity (",Slab[2],")"),
          zlab=c("Characteristic Life (",Llab,")"))
    #fig<-plot_ly(x = ~S1line, y = ~S2line, z = ~life(MLEoutput[[1]],Sline), type = 'mesh3d')
    #fig <- fig %>% layout(scene = list(xaxis = list(title = c("Characteristic Stress (",Slab[1],")")),
    #                                   yaxis = list(title = c("Characteristic Humidity (",Slab[2],")")),
    #                                   zaxis = list(title = c("Characteristic Life (",Llab,")"))))
    #fig
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lifeUse <- function(theta,SF) {
      theta[ishift+3]/((SF[2]^theta[ishift+2])*exp(-theta[ishift+1]/SF[1]))
    } # Single stress entry life (for Use Life primarily)
    life <- function(theta,SF1,SF2) {
      theta[ishift+3]/((SF2^theta[ishift+2])*exp(-theta[ishift+1]/SF1))
    } # Multiple stress entry life (for vectors)
    UselifeLSQ <- lifeUse(LSQoutput[[3]],SUse)
    UselifeMLE <- lifeUse(MLEoutput[[1]],SUse)

  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    life <- function(theta,SF) {
      (1/SF[1])*exp((theta[ishift+1] + (theta[ishift+2]/SF[1])) + (theta[ishift+3] + (theta[ishift+4]/SF[1]))*SF[2])
    }

  }

  # Compute Use life
  #UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #UselifeMLE <- life(MLEoutput[[1]],SUse)

  # Return parameter list
  return(list(UselifeLSQ,UselifeMLE))
}
