# Fail and Censor Count
# Developed by Dr. Reuel Smith, 2020-2022

count.failcen <- function(mat,n) {
  # Removes all duplicates from xi column
  xu <- mat[,1][!duplicated(mat[,1])]
  # Initialize the output matrix where column 2 is the number of
  # failures at  x_i, column 3 is the number of censored units
  # at x_i, and column 4 is the remaining units at x_i
  #
  matfailcencount<-matrix(c(xu,rep(0, length(xu)),rep(0, length(xu)),n,rep(0, length(xu)-1)), nrow = length(xu), ncol = 4)
  i3<-1
  for (i2 in 1:n) {
    if (mat[i2,1]==xu[i3]) {
      i3<-i3
    } else {
      i3<-i3+1
      matfailcencount[i3,4]<-n-sum(matfailcencount[1:i3-1,2])-sum(matfailcencount[1:i3-1,3])
    }
    if (mat[i2,2]==1 && i3<=n) {
      matfailcencount[i3,2]<-matfailcencount[i3,2]+1
    }
    if (mat[i2,2]==0 && i3<=n) {
      matfailcencount[i3,3]<-matfailcencount[i3,3]+1
    }
  }
  return(matfailcencount)
}
