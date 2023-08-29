# Life from Life Strain Properties
# Developed by Reuel Smith, 2022

lifestrain.life.trace <- function(strain_amp,E,b,c,sig_f,eps_f,stressunits,options){
  # dat is entered as a list made up of stress, strain, and cycles in that order
  library(pracma)

  if(missing(E)){
    stop('Enter modulus of elasticity E (in MPa or ksi).')
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT" && length(options$stress_amp) == 0){
    stop('SWT can only be computed with current load condition of stress amplitude (stress_amp) entered in options.')
  }

  # Check units and set up axis labels
  if(missing(stressunits) || stressunits == 1){
    stressunits <- c("MPa")
  }
  if(stressunits == 2){
    stressunits <- c("ksi")
  }

  # Find Cycles and Reversals to Failure
  if(missing(options)==TRUE || (isFALSE(missing(options)) && length(options$mean_stress_corr)==0)){
    # Coffin-Mason
    straintrace <- function(reversal) strain_amp - (sig_f/E)*((reversal)^b) - eps_f*((reversal)^c)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Morrow" ){
    # Morrow Mean Stress Correction
    Sm <- options$Sm
    straintrace <- function(reversal) strain_amp - ((sig_f - Sm)/E)*((reversal)^b) - eps_f*((reversal)^c)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "ModifiedMorrow" ){
    # Modified Morrow Mean Stress Correction
    Sm <- options$Sm
    straintrace <- function(reversal) strain_amp - ((sig_f - Sm)/E)*((reversal)^b) - eps_f*(((sig_f - Sm)/sig_f)^(c/b))*((reversal)^c)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT" && length(options$stress_amp) == 1){
    # SWT Mean Stress Correction
    Sm <- options$Sm
    stress_amp <- options$stress_amp
    stressmax_corr <- stress_amp + Sm
    straintrace <- function(reversal) strain_amp - (1/stressmax_corr)*((sig_f^2)/E)*((reversal)^(2*b)) - (1/stressmax_corr)*eps_f*sig_f*((reversal)^(b+c))
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Walker"  && (length(options$R) == 1 || length(options$stressrange) == 2) &&  length(options$gam) == 1){
    # Walker Mean Stress Correction
    if(length(options$R) == 0 && length(options$stressrange) == 2){
      R <- min(options$stressrange)/max(options$stressrange)
    }
    if(length(options$R) == 1 && length(options$stressrange) == 0){
      R <- options$R
    }
    gam <- options$gam

    straintrace <- function(reversal) strain_amp - (sig_f/E)*((0.5*(1 - R))^(1 - gam))*((reversal)^b) - eps_f*((0.5*(1 - R))^((c/b)*(1 - gam)))*((reversal)^c)
  }

  rev_2N_trace <- findzeros(straintrace, 10, 10^9)
  Life_trace <- rev_2N_trace/2

  return(list(strainreversals = rev_2N_trace,strainlife = Life_trace))
}
