# Tukey Non-Parametric Output Tabulation
# Developed by Dr. Reuel Smith, 2020-2022

plotposit.tukey <- function(i, xi, rc) {
  # Obtain full length of all data and cut data
  nx<-nxi.count(xi,rc)
  # Get differentials of xi
  dx <- c(diff(nx[[2]]),0)
  # Tukey plot for F(x) and R(x)
  F<-(i - (1/3))/(nx[[1]]+(1/3))
  R<-1-F
  H<- -log(R)
  h<-1/((nx[[1]]-i+(2/3))*dx)
  f<-1/((nx[[1]]+(1/3))*dx)
  matFR<-matrix(c(nx[[2]],F,R,h,H,f), nrow = length(nx[[2]]), ncol = 6)
  return(matFR)
}
