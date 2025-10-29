# Failure Time/Right-Censored Data sort (Single-Stress)
# Developed by Dr. Reuel Smith, 2020-2022

sort.xircdata <- function(rawdat) {
  # Input has to be a two column table
  if(length(rawdat[,2])==sum(rawdat[,2])) { # No right censored data case
    xi <- sort(rawdat[,1]) # sort the xi (event) data from least to greatest
    i <- rankcalc(xi)  # Rank data with `rankcalc` function
    xitab<-list(i, xi, NULL) # Retabulate data
  }
  else { # right censored data case
    xi <- sort(rawdat[,1][which(rawdat[,2] == 1)]) # sort the xi (event) data from least to greatest
    rc <- sort(rawdat[,1][which(rawdat[,2] == 0)]) # sort the rc (censored event) data from least to greatest
    i <- rankcalc(xi,rc) # Rank data with `rankcalc` function
    xitab<-list(i,xi,rc) # Retabulate data
  }
  return(xitab) # return retabulated data
}
