# Bernard Bos-Levenbach Non-Parametric Output Tabulation
# Developed by Dr. Reuel Smith, 2020-2022

plotposit.bernbos <- function(i, xi, rc) {
  # Obtain full length of all data and cut data
  nx<-nxi.count(xi,rc)
  # Bernard Bos-Levenbach plot for F(x) and R(x)
  F<-(i - 0.3)/(nx[[1]]+0.2)
  R<-1-F
  H<- -log(R)
  h<-c(H[1],diff(H))
  matFR<-matrix(c(nx[[2]],F,R,h,H), nrow = length(nx[[2]]), ncol = 5)
  return(matFR)
}
