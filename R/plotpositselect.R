# Plotting Position Selector
# Developed by Dr. Reuel Smith, 2020-2022

plotposit.select <- function(xi, rc, pp) {
  if (pp=="Blom") {
    if(is.null(rc)) {
      matFR<-plotposit.blom(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.blom(rankcalc(xi,rc),xi,rc)
      }
  }
  if (pp=="Mean") {
    if(is.null(rc)) {
      matFR<-plotposit.mean(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.mean(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="Median") {
    if(is.null(rc)) {
      matFR<-plotposit.median(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.median(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="Midpoint") {
    if(is.null(rc)) {
      matFR<-plotposit.midpt(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.midpt(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="Beard") {
    if(is.null(rc)) {
      matFR<-plotposit.beard(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.beard(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="BernardBosLevenbach") {
    if(is.null(rc)) {
      matFR<-plotposit.bernbos(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.bernbos(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="Tukey") {
    if(is.null(rc)) {
      matFR<-plotposit.tukey(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.tukey(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="Grigorten") {
    if(is.null(rc)) {
      matFR<-plotposit.grigorten(rankcalc(xi),xi)
    }
    else {
      matFR<-plotposit.grigorten(rankcalc(xi,rc),xi,rc)
    }
  }
  if (pp=="NelsonAalen") {
    if(is.null(rc)) {
      matFR<-plotposit.nelsonaalen(xi)
    }
    else {
      matFR<-plotposit.nelsonaalen(xi,rc)
    }
  }
  if (pp=="KaplanMeier") {
    if(is.null(rc)) {
      matFR<-plotposit.kaplanmeier(xi)
    }
    else {
      matFR<-plotposit.kaplanmeier(xi,rc)
    }
  }
  return(matFR)
}
