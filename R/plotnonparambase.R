# Non-Parametric Plot Fitter
# Developed by Dr. Reuel Smith, 2020-2022

plot.nonparam <- function(xi, yi) {
  n<-length(xi)
  n2<-2*(1+n)
  xout<-rep(0, n2)

  if(yi[1]<yi[n] | yi[1]<yi[n-1]) {
    # Plot regeneration for unreliability data
    yout<-rep(0, n2)
  } else {
    # Plot regeneration for reliability data
    yout<-rep(1, n2)
  }
  # build x and y output vector
  for(i in 1:n) {
    xout[2*i]<-xi[i]
    xout[2*i+1]<-xi[i]
    yout[2*i+1]<-yi[i]
    yout[2*i+2]<-yi[i]
  }
  xout[n2]<-(round(xi[n]/10^round(log10(xi[n])))+1)*10^round(log10(xi[n]))
  return(list(xout,yout))
}
