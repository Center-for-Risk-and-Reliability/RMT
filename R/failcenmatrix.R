# Fail and Censor Matrix
# Developed by Dr. Reuel Smith, 2020-2022

matrix.failcen <- function(xi, rc, nx) {
  # Create ones vector for uncensored units
  failed<-rep(1, length(xi))
  # Sort and order failure and censored data
  if(missing(rc)) {
    faicen<-matrix(c(xi,failed), nrow = nx, ncol = 2)
  } else {
    # Create zeroes vector for censored units
    cen<-rep(0, length(rc))
    faicen<-matrix(c(xi,rc,failed,cen), nrow = nx, ncol = 2)
  }
  faicenmat<-faicen[order(faicen[,1]),]
  return(faicenmat)
  }
