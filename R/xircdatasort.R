# Failure Time/Right-Censored Data sort (Single-Stress)
# Developed by Dr. Reuel Smith, 2020-2022

sort.xircdata <- function(rawdat) {
  # Input has to be a two column table
  if(length(rawdat[,2])==sum(rawdat[,2])) {
    xi <- sort(rawdat[,1])
    i <- rankcalc(xi)
    xitab<-list(i, xi, NULL)
  }
  else {
    xi <- sort(rawdat[,1][which(rawdat[,2] == 1)])
    rc <- sort(rawdat[,1][which(rawdat[,2] == 0)])
    i <- rankcalc(xi,rc)
    xitab<-list(i,xi,rc)
  }
  return(xitab)
}
