# Kaplan-Meier Non-Parametric Output Tabulation
# Developed by Dr. Reuel Smith, 2020-2022

plotposit.kaplanmeier <- function(xi, rc) {
  # Obtain full length of all data and cut data
  nx<-nxi.count(xi,rc)
  # Form Matrix
  if (missing(rc)){
    faicen<-matrix.failcen(xi, nx=nx[[1]])
  } else {
    faicen<-matrix.failcen(xi, rc, nx=nx[[1]])
  }
  # Get the number of failures and censored units
  fail_and_cen<-count.failcen(faicen,nx[[1]])
  # Isolate the number of failures
  di<-fail_and_cen[,2][which(fail_and_cen[,2] > 0)]
  # Isolate the number of units remaining
  ni<-fail_and_cen[,4][which(fail_and_cen[,2] > 0)]
  # Kaplan-Meier plot for F(x) and R(x)
  R <- cumprod((ni-di)/ni)
  F <- 1-R
  H <- -log(R)
  h <- -log((ni-di)/ni)
  f <- h*R
  matFR<-matrix(c(nx[[2]],F,R,h,H,f), nrow = length(nx[[2]]), ncol = 6)
  return(matFR)
}
