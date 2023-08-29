# Non-Parametric Output Tabulation
# Developed by Dr. Reuel Smith, 2020-2022

plottable.nonparam <- function(xi, rc, FRhH, relfcn, alpha, xlabel) {
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
  bounds<-conf.greenwood(alpha,nxi,FRhH[,3],count.failcen(rawdat,nxi))
  if(relfcn=="unreliability"){
    medplot<-plot.nonparam(FRhH[,1], FRhH[,2])
    lowplot<-plot.nonparam(FRhH[,1], bounds[[3]])
    hiplot<-plot.nonparam(FRhH[,1], bounds[[4]])
    plot(medplot[[1]], medplot[[2]], type="l", xlab=xlabel, ylab="Unreliability",ylim=c(0,1) , col="blue")
    lines(medplot[[1]],lowplot[[2]], type="l", lty=2, col="blue")
    lines(medplot[[1]],hiplot[[2]], type="l", lty=2, col="blue")
  }
  if(relfcn=="reliability"){
    medplot<-plot.nonparam(FRhH[,1], FRhH[,3])
    lowplot<-plot.nonparam(FRhH[,1], bounds[[1]])
    hiplot<-plot.nonparam(FRhH[,1], bounds[[2]])
    plot(medplot[[1]], medplot[[2]], type="l", xlab=xlabel, ylab="Reliability",ylim=c(0,1) , col="blue")
    lines(medplot[[1]],lowplot[[2]], type="l", lty=2, col="blue")
    lines(medplot[[1]],hiplot[[2]], type="l", lty=2, col="blue")
  }
  if(relfcn=="hazard"){
    medplot<-plot.nonparamhaz(FRhH[,1], FRhH[,4])
    lowplot<-plot.nonparamhaz(FRhH[,1], bounds[[7]])
    hiplot<-plot.nonparamhaz(FRhH[,1], bounds[[8]])
    if (is.finite(max(hiplot[[2]]))){
      ylims<-c(0,max(hiplot[[2]]))
    } else {
      ylims<-c(0,hiplot[[2]][length(hiplot[[2]])-2])
    }
    plot(medplot[[1]], medplot[[2]], type="l", xlab=xlabel, ylab="Hazard",ylim=ylims, col="blue")
    lines(medplot[[1]],lowplot[[2]], type="l", lty=2, col="blue")
    lines(medplot[[1]],hiplot[[2]], type="l", lty=2, col="blue")
  }
  if(relfcn=="cumulativehazard"){
    medplot<-plot.nonparamhaz(FRhH[,1], FRhH[,5])
    lowplot<-plot.nonparamhaz(FRhH[,1], bounds[[5]])
    hiplot<-plot.nonparamhaz(FRhH[,1], bounds[[6]])
    if (is.finite(max(hiplot[[2]]))){
      ylims<-c(0,max(hiplot[[2]]))
    } else {
      ylims<-c(0,hiplot[[2]][length(hiplot[[2]])-2])
    }
    plot(medplot[[1]], medplot[[2]], type="l", xlab=xlabel, ylab="Cumulative Hazard",ylim=ylims, col="blue")
    lines(medplot[[1]],lowplot[[2]], type="l", lty=2, col="blue")
    lines(medplot[[1]],hiplot[[2]], type="l", lty=2, col="blue")
  }
  nonparamset<-matrix(c(FRhH[,1],FRhH[,2],bounds[[3]],bounds[[4]],FRhH[,3],bounds[[1]],bounds[[2]],FRhH[,4],bounds[[7]],bounds[[8]],FRhH[,5],bounds[[5]],bounds[[6]]),nrow = length(xi), ncol = 13)
  colnames(nonparamset)<-c(xlabel,"Unreliability","Low bound","High bound","Reliability","Low bound","High bound","Hazard","Low bound","High bound","Cumulative Hazard","Low bound","High bound")
  return(nonparamset)
}
