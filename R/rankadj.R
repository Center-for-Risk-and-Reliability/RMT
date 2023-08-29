# Rank Sort by Adjusted Rank Calculation
# Developed by Dr. Reuel Smith, 2021-2022

rankadj <- function(xi, rc) {
  # Obtain full length of all data
  if(missing(rc)) {
    nxitot<-nxi.count(xi)
  } else {
    nxitot<-nxi.count(xi,rc)
  }
  # Form Matrix
  faicen<-matrix.failcen(xi, rc, nx=nxitot[[1]])
  # Get the number of failures and censored units
  fail_and_cen<-count.failcen(faicen,nxitot[[1]])
  # Isolate the number of failures
  f<-fail_and_cen[,2][which(fail_and_cen[,2] > 0)]
  # Isolate the number of censored
  c<-fail_and_cen[,3][which(fail_and_cen[,2] > 0)]
  # Compute duplicates m
  m<-f+c-1
  # Get indicies of failure
  og_i<-which(duplicated(faicen[,1]) %in% FALSE)[which(fail_and_cen[,2] > 0)]
  iadj<-og_i + 0.5*(m-c)
  # Calculate rank adjusted rank
  i<-rep(1, length(iadj))
  if(identical(iadj,1) && length(iadj) == 1){
    i <- 1
  } else {
    # Only applies for iadj >= 2
    for (i2 in 2:length(iadj)) {
      i[i2]<-i[i2 - 1] + ((nxitot[[1]] + 1)-i[i2 - 1])/(2 + nxitot[[1]] - og_i[i2])
    }
  }

  return(i)
}
