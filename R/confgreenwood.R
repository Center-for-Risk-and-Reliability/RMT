# Non-Parametric Greenwood Confidence
# Developed by Dr. Reuel Smith, 2020-2022

conf.greenwood <- function(alp,n,R,nandd) {
  # Trim the failures (di) and remaining units (ni) of only needed cases
  di<-nandd[,2][which(nandd[,2] > 0)]
  ni<-nandd[,4][which(nandd[,2] > 0)]
  # Compute
  Rlow<-R-qt(1-alp/2, df=n-1)*sqrt((R^2)*cumsum(di/(ni*(ni-di))))
  Rhi<-R+qt(1-alp/2, df=n-1)*sqrt((R^2)*cumsum(di/(ni*(ni-di))))
  Hlow<- -log(R)-qt(1-alp/2, df=n-1)*sqrt(cumsum(di/(ni*(ni-di))))
  Hhi<- -log(R)+qt(1-alp/2, df=n-1)*sqrt(cumsum(di/(ni*(ni-di))))
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
