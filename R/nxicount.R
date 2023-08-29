# Data Count by Fail Time and Right-Censored
# Developed by Dr. Reuel Smith, 2021-2022

nxi.count <- function(xi, rc) {
  # Check if there is right censored data or not
  if(missing(rc)) {
    # Obtain full length of all data
    n<-sum(length(xi))
  }
  else {
    # Obtain full length of all data
    n<-sum(length(xi),length(rc))
  }
  # Cut and sort xi
  xired <- sort(xi[!duplicated(xi)])
  return(list(n,xired))
}
