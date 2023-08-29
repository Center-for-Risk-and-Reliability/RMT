# Calculate Quantiles by Column
# Developed by Dr. Reuel Smith, 2021-2022

colQuantiles <- function(mat,prob){
  # NOTE: Only works for one probability right now.  Will work on multiple later.
  q<-rep(0,dim(mat)[2])
  for (i in 1:dim(mat)[2]) {
    q[i]<-as.vector(quantile(mat[,i],prob))
  }
  return(q)
}
