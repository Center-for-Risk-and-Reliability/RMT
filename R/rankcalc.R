# Rank Calculator
# Developed by Dr. Reuel Smith, 2021-2022

rankcalc <- function(xi, rc) {
  # Sort failure data xi
  xi<-sort(xi);
  # Check if there is right censored data or not
  if(missing(rc)) {
    # Compute the rank for a non-censored case
    i <- rank(xi)
    # Remove duplicate ranks for duplicate xi
    i <- i[!duplicated(i)]
    }
  else {
    # run rank adjustment for else case
    i <- rankadj(xi, rc)
  }
  return(i)
}
