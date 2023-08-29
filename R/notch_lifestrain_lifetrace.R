# Life from Life Strain Properties with Notch Effect
# Developed by Reuel Smith, 2022

notch.lifestrain.life.trace <- function(dimensions,geometry,E,b,c,sig_f,eps_f,K,n,stressunits = 1,loadconditions,options){
  library(pracma)

  if(missing(E)){
    stop('Enter modulus of elasticity E (in MPa or ksi).')
  }
  # if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT" && length(options$stress_amp) == 0){
  #   stop('SWT can only be computed with current load condition of stress amplitude (stress_amp) entered in options.')
  # }

  # Check units and set up axis labels
  if(stressunits == 1){
    stressunitslabel <- c("MPa")
  }
  if(stressunits == 2){
    stressunitslabel <- c("ksi")
  }

  # Compute the stress concentration factor, Kt
  Kt <- stress.concentration.factor(dimensions,geometry)

  # Calculate Life to Failure based on sample type and area
  if(geometry == "rect_1semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$t) == 1 && length(dimensions$r) == 1){
    Area <- dimensions$t*(dimensions$W - dimensions$r)
  }
  if(geometry == "rect_2semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$t) == 1 && length(dimensions$r) == 1){
    Area <- dimensions$t*(dimensions$W - 2*dimensions$r)
  }
  Snet <- loadconditions$Sa/Area

  if(Kt*Snet > loadconditions$Sy){
    # Compute the notch sensitivity factor, q
    # metric (mm)
    if(length(loadconditions$a)!=1 && stressunits == 1){
      a <- 0.0254*(((300*6.8947572932)/loadconditions$Su)^1.8)
    }
    # English (inches)
    if(length(loadconditions$a)==0 && stressunits == 2){
      a <- (1e-3)*((300/loadconditions$Su)^1.8)
    }
    if(length(loadconditions$a)==1){
      a <- loadconditions$a
    }
    q <- 1/(1 + (a/dimensions$r))

    # Compute fatigue stress concentration factor, Kf
    Kf <- 1 + q*(Kt - 1)

    # Compute maxstress and maxstrain value
    stress_strain <- ((Kf*Snet)^2)/E
    stresstrace <- function(maxstress) (maxstress/E) + ((maxstress/K)^(1/n)) - stress_strain/maxstress
    stressmax <- findzeros(stresstrace, 10, E)
    strainmax <- stress_strain/stressmax

    # Compute minstress and minstrain value
    Dstress_Dstrain <- ((Kf*2*Snet)^2)/E
    Dstresstrace <- function(Dstress) (Dstress/E) + 2*((Dstress/(2*K))^(1/n)) - Dstress_Dstrain/Dstress
    Dstress <- findzeros(Dstresstrace, 10, E)
    Dstrain <- Dstress_Dstrain/Dstress
    stressmin <- stressmax - Dstress
    strainmin <- strainmax - Dstrain
  }
  Sm <- 0.5*(stressmax + stressmin)
  Stress_amp <- 0.5*(stressmax - stressmin)
  strain_amp <- 0.5*(strainmax - strainmin)



  # Find Cycles and Reversals to Failure
  if(missing(options)==TRUE || (isFALSE(missing(options)) && length(options$mean_stress_corr)==0)){
    # Coffin-Mason
    straintrace <- function(reversal) strain_amp - (sig_f/E)*((reversal)^b) - eps_f*((reversal)^c)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Morrow" ){
    # Morrow Mean Stress Correction
    straintrace <- function(reversal) strain_amp - ((sig_f - Sm)/E)*((reversal)^b) - eps_f*((reversal)^c)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "ModifiedMorrow" ){
    # Modified Morrow Mean Stress Correction
    straintrace <- function(reversal) strain_amp - ((sig_f - Sm)/E)*((reversal)^b) - eps_f*(((sig_f - Sm)/sig_f)^(c/b))*((reversal)^c)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT"){
    # SWT Mean Stress Correction
    stressmax_corr <- stress_amp + Sm
    straintrace <- function(reversal) strain_amp - (1/stressmax_corr)*((sig_f^2)/E)*((reversal)^(2*b)) - (1/stressmax_corr)*eps_f*sig_f*((reversal)^(b+c))
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Walker"  &&  length(options$gam) == 1){
    # Walker Mean Stress Correction
    gam <- options$gam

    straintrace <- function(reversal) strain_amp - (sig_f/E)*((0.5*(1 - R))^(1 - gam))*((reversal)^b) - eps_f*((0.5*(1 - R))^((c/b)*(1 - gam)))*((reversal)^c)
  }

  rev_2N_trace <- findzeros(straintrace, 10, 10^9)
  Life_trace <- rev_2N_trace/2

  return(list(Kf = Kf, stress_strain = stress_strain, stressmax = stressmax, strainmax = strainmax, stressmin = stressmin, strainmin = strainmin, strainreversals = rev_2N_trace,strainlife = Life_trace))
}
