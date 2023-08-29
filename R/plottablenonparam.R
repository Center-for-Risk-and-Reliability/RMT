# Non-Parametric Output Tabulation
# Developed by Dr. Reuel Smith, 2020-2022

plottable.nonparam <- function(xi, rc, FRhH, relfcn, alpha, xlabel) {
  library(ggplot2)
  if (missing(rc)){
    nxi<-length(xi)
    rawdat<-matrix.failcen(xi, nx=nxi)
  } else {
    nxi<-length(c(xi,rc))
    rawdat<-matrix.failcen(xi, rc, nx=nxi)
  }
  if (missing(xlabel)){
    xlabel<-"X"
  }

  bounds<-conf.binomial(alpha,nxi,FRhH[,3])
  if(relfcn=="probabilitydensity"){
    X <- FRhH[,1]
    dens_X <- FRhH[,6][1:(length(FRhH[,1])-1)]
    df <- data.frame(X = FRhH[,1][1:(length(FRhH[,1])-1)], Fx = FRhH[,6][1:(length(FRhH[,1])-1)])
    plotout<-ggplot() +
      geom_col(data=df, aes(X,Fx), colour = 'blue', size = 1.9) +
      xlab(xlabel) +
      ylab("Probability density")
  }
  if(relfcn=="unreliability"){
    medplot<-plot.nonparam(FRhH[,1], FRhH[,2])
    lowplot<-plot.nonparam(FRhH[,1], bounds[[3]])
    hiplot<-plot.nonparam(FRhH[,1], bounds[[4]])
    df <- data.frame(X = medplot[[1]], Fx = medplot[[2]], Fx_low = lowplot[[2]], Fx_hi = hiplot[[2]])
    plotout<-ggplot() +
      geom_line(data=df, aes(X,Fx), colour = 'blue', size = 0.9, linetype = "dashed") +
      xlab(xlabel) +
      ylab("Unreliability")
    plotout <- plotout + geom_ribbon(data=df, aes(ymin=Fx_low, ymax=Fx_hi, x=X), alpha=0.5, fill = "red")
  }
  if(relfcn=="reliability"){
    medplot<-plot.nonparam(FRhH[,1], FRhH[,3])
    lowplot<-plot.nonparam(FRhH[,1], bounds[[1]])
    hiplot<-plot.nonparam(FRhH[,1], bounds[[2]])
    df <- data.frame(X = medplot[[1]], Fx = medplot[[2]], Fx_low = lowplot[[2]], Fx_hi = hiplot[[2]])
    plotout<-ggplot() +
      geom_line(data=df, aes(X,Fx), colour = 'blue', size = 0.9, linetype = "dashed") +
      xlab(xlabel) +
      ylab("Reliability")
    plotout <- plotout + geom_ribbon(data=df, aes(ymin=Fx_low, ymax=Fx_hi, x=X), alpha=0.5, fill = "red")
  }
  if(relfcn=="hazard"){
    X <- FRhH[,1]
    dens_X <- FRhH[,4][1:(length(FRhH[,1])-1)]
    df <- data.frame(X = FRhH[,1][1:(length(FRhH[,1])-1)], Fx = FRhH[,4][1:(length(FRhH[,1])-1)])
    plotout<-ggplot() +
      geom_col(data=df, aes(X,Fx), colour = 'red', size = 1.9) +
      xlab(xlabel) +
      ylab("Hazard")
  }
  if(relfcn=="cumulativehazard"){
    medplot<-plot.nonparamhaz(FRhH[,1], FRhH[,5])
    lowplot<-plot.nonparamhaz(FRhH[,1], bounds[[5]])
    hiplot<-plot.nonparamhaz(FRhH[,1], bounds[[6]])
    df <- data.frame(X = medplot[[1]], Fx = medplot[[2]], Fx_low = lowplot[[2]], Fx_hi = hiplot[[2]])
    plotout<-ggplot() +
      geom_line(data=df, aes(X,Fx), colour = 'blue', size = 0.9, linetype = "dashed") +
      xlab(xlabel) +
      ylab("Cumulative Hazard")
    plotout <- plotout + geom_ribbon(data=df, aes(ymin=Fx_low, ymax=Fx_hi, x=X), alpha=0.5, fill = "red")
  }
  nonparamset<-matrix(c(FRhH[,1],FRhH[,6],FRhH[,2],bounds[[3]],bounds[[4]],FRhH[,3],bounds[[1]],bounds[[2]],FRhH[,4],FRhH[,5],bounds[[5]],bounds[[6]]),nrow = length(xi), ncol = 12)
  colnames(nonparamset)<-c(xlabel,"Probability Density","Unreliability","Low bound","High bound","Reliability","Low bound","High bound","Hazard","Cumulative Hazard","Low bound","High bound")
  return(list(nonparamplot = plotout, paramtable = nonparamset))
}
