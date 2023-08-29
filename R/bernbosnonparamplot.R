# Bernard Bos-Levenbach Non-Parametric Plot
# Developed by Dr. Reuel Smith, 2020-2022

nonparamplot.bernbos <- function(xi, rc, relfcn, confbnd, xlabel) {
  if(missing(relfcn)) {
    relfcn<-"reliability"
  }

  if(missing(confbnd)) {
    alpha<-1-0.95
  } else {
    alpha<-(100-confbnd)*0.01
  }
  if(missing(xlabel)) {
    xlabel<-"X"
  }

  # Compute the unreliability, reliability, hazard, and cumulative hazard

  if (missing(rc)){
    FRhH<-plotposit.bernbos(rankcalc(xi),xi)
    nonparamset<-plottable.nonparam(xi, NULL, FRhH, relfcn, alpha, xlabel)
  } else {
    FRhH<-plotposit.bernbos(rankcalc(xi,rc),xi, rc)
    nonparamset<-plottable.nonparam(xi, rc, FRhH, relfcn, alpha, xlabel)
  }


  return(nonparamset)
}
