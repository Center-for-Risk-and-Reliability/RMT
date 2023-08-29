# Non-Parametric Binomial Confidence
# Developed by Dr. Reuel Smith, 2020-2022

conf.binomial <- function(alp,n,R) {

  # Compute
  Fmean <- 1 - R
  confid <- 1 - alp
  Z <- qnorm(1-(1-confid)/2,0,1)
  Fdiff <- Z*sqrt((R*(1-R))/n)
  Fhi <- Fmean + Fdiff
  Fhi[which(Fhi>0.9999)] <- 0.9999
  Flow <- Fmean - Fdiff
  Flow[which(Flow<0.0001)] <- 0.0001

  Rlow<-1-Fhi
  Rhi<-1-Flow

  Hlow<- -log(Rhi)
  Hhi<- -log(Rlow)
  # Hlow<- -log(R) + log(1 - Fdiff)
  # Hhi<- -log(R) - log(1 - Fdiff)

  # Hlow<- -log(R)-qt(1-alp/2, df=n-1)*sqrt((1-R)/R)
  # Hhi<- -log(R)+qt(1-alp/2, df=n-1)*sqrt((1-R)/R)


  Hlow[which(Hlow<0)]<-0

  hlow<-c(Hlow[1],diff(Hlow))
  hlow[which(hlow<0)]<-0
  hhi<-c(Hhi[1],diff(Hhi))

  # Group
  return(list(Rlow,Rhi,Flow,Fhi,Hlow,Hhi,hlow,hhi))
}
