# Plotting Position Selector
# Developed by Dr. Reuel Smith, 2020-2025

plotposit.select <- function(xi, rc=NULL, pp="Blom",CDFgiven = NULL) {
  # RCS (2/25/2025) - Added CDFgiven vector input to compute rank given xi and the percent failure
  if (pp=="Blom") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.blom(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.blom(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + 0.25) + 0.375
      matFR<-plotposit.blom(rank_new,xi)
    }
  }
  if (pp=="Mean") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.mean(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.mean(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + 1)
      matFR<-plotposit.mean(rank_new,xi)
    }
  }
  if (pp=="Median") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.median(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.median(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + 0.4) + 0.3
      matFR<-plotposit.median(rank_new,xi)
    }
  }
  if (pp=="Midpoint") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.midpt(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.midpt(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*length(xi) + 0.5
      matFR<-plotposit.midpt(rank_new,xi)
    }
  }
  if (pp=="Beard") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.beard(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.beard(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + 0.38) + 0.31
      matFR<-plotposit.beard(rank_new,xi)
    }
  }
  if (pp=="BernardBosLevenbach") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.bernbos(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.bernbos(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + 0.2) + 0.3
      matFR<-plotposit.bernbos(rank_new,xi)
    }
  }
  if (pp=="Tukey") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.tukey(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.tukey(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + (1/3)) + (1/3)
      matFR<-plotposit.tukey(rank_new,xi)
    }
  }
  if (pp=="Grigorten") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.grigorten(rankcalc(xi),xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.grigorten(rankcalc(xi,rc),xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      rank_new <- CDFgiven*(length(xi) + 0.12) + 0.44
      matFR<-plotposit.grigorten(rank_new,xi)
    }
  }
  if (pp=="NelsonAalen") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.nelsonaalen(xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.nelsonaalen(xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      matFR<-plotposit.nelsonaalen(xi)
      matFR[,2] <- CDFgiven
      matFR[,3] <- 1-CDFgiven
      matFR[,5] <- -log(1-CDFgiven)
    }
  }
  if (pp=="KaplanMeier") {
    if(is.null(rc) == TRUE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.kaplanmeier(xi)
    }
    if(is.null(rc) == FALSE && is.null(CDFgiven) == TRUE) {
      matFR<-plotposit.kaplanmeier(xi,rc)
    }
    if(is.null(CDFgiven) == FALSE) {
      matFR<-plotposit.kaplanmeier(xi)
      matFR[,2] <- CDFgiven
      matFR[,3] <- 1-CDFgiven
      matFR[,5] <- -log(1-CDFgiven)
    }
  }
  return(matFR)
}
