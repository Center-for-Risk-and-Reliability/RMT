# Nelson-Aalen Non-Parametric Output Tabulation
# Developed by Dr. Reuel Smith, 2020-2022

plotposit.nelsonaalen <- function(xi, rc) {
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
  # Nelson-Aalen plot for F(x) and R(x)
  F<-1-exp(-cumsum(di/ni))
  R<-1-F
  H<- -log(R)
  h<-c(H[1],diff(H))
  matFR<-matrix(c(nx[[2]],F,R,h,H), nrow = length(nx[[2]]), ncol = 5)
  return(matFR)
}
