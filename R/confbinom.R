# Non-Parametric Binomial Confidence
# Developed by Dr. Reuel Smith, 2020-2022

conf.binomial <- function(alp,n,R) {

  # Compute
  Rlow<-R-qt(1-alp/2, df=n-1)*sqrt(R*(1-R))
  Rhi<-R+qt(1-alp/2, df=n-1)*sqrt(R*(1-R))
  Hlow<- -log(R)-qt(1-alp/2, df=n-1)*sqrt((1-R)/R)
  Hhi<- -log(R)+qt(1-alp/2, df=n-1)*sqrt((1-R)/R)
  # Rounds highs and lows greater than 1 and less than 0 respectively
  Rlow[which(Rlow<0)]<-0
  Rlow[which(is.nan(Rlow))]<-0
  Rhi[which(Rhi>1)]<-1
  Rhi[which(is.nan(Rhi))]<-0


  Flow<-1-Rhi
  Fhi<-1-Rlow

  Hlow[which(Hlow<0)]<-0

  hlow<-c(Hlow[1],diff(Hlow))
  hlow[which(hlow<0)]<-0
  hhi<-c(Hhi[1],diff(Hhi))

  # Group
  return(list(Rlow,Rhi,Flow,Fhi,Hlow,Hhi,hlow,hhi))
}
