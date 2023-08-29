# Non-Parametric Hazard Plot Fitter
# Developed by Dr. Reuel Smith, 2020-2022

plot.nonparamhaz <- function(xi, yi) {
  n<-length(xi)
  n2<-2*(1+n)
  xout<-rep(0, n2)
  yout<-rep(0, n2)
  # build x and y output vector
  for(i in 1:n) {
    xout[2*i]<-xi[i]
    xout[2*i+1]<-xi[i]
    yout[2*i+1]<-yi[i]
    yout[2*i+2]<-yi[i]
  }
  xout[n2]<-(round(xi[n]/10)+1)*10
  return(list(xout,yout))
}
