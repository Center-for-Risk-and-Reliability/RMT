# Life from Life Strain Properties with Notch Effect
# Developed by Reuel Smith, 2022

notch.lifestrain.life.trace <- function(dimensions,geometry=list(r = NULL, R = NULL, W = NULL, d = NULL, D = NULL, alp = NULL, a = NULL, b = NULL, W_m_r = NULL, W_m_d = NULL, W_m_2d = NULL, W_m_2r = NULL, t = NULL),E,b,c,sig_f,eps_f,K,n,stressunits = 1,loadconditions = list(Su = NULL, Sy = NULL, Sa = NULL, Sa1 = NULL, Sa2 = NULL, p = NULL, a = NULL,  q = NULL),options = list(mean_stress_corr = NULL, Kf = NULL)){
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

  # Calculate Life to Failure based on sample type and area
  if(geometry == "rect_1semicirc_edge" || (geometry == "rect_1U_edge" && dimensions$r == dimensions$d)){
    if(is.null(dimensions$W)==FALSE && is.null(dimensions$r)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*(dimensions$W - dimensions$r)
    }
    if(is.null(dimensions$W_m_r)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*dimensions$W_m_r
    }
  }
  if(geometry == "rect_2semicirc_edge" || geometry == "plate_circ_offset_hole" || geometry == "rect_inf_semicirc_edge" || geometry == "rect_double_fillet"){
    if(is.null(dimensions$W)==FALSE && is.null(dimensions$r)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*(dimensions$W - 2*dimensions$r)
    }
    if(is.null(dimensions$W_m_2r)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*dimensions$W_m_2r
    }
  }
  if((geometry == "rect_1U_edge" && dimensions$r != dimensions$d) || geometry == "rect_circ_symm_hole" || geometry == "rect_ellips_symm_hole"){
    if(is.null(dimensions$W)==FALSE && is.null(dimensions$d)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*(dimensions$W - dimensions$d)
    }
    if(is.null(dimensions$W_m_d)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*dimensions$W_m_d
    }
  }
  if((geometry == "rect_2U_edge" && dimensions$r != dimensions$d) || geometry == "rect_2V_edge"){
    if(is.null(dimensions$W)==FALSE && is.null(dimensions$r)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*(dimensions$W - 2*dimensions$d)
    }
    if(is.null(dimensions$W_m_2r)==FALSE && is.null(dimensions$t)==FALSE){
      Area <- dimensions$t*dimensions$W_m_2d
    }
  }
  if(geometry == "shaft_fillet_single" || geometry == "shaft_reg_U_groove" || geometry == "shaft_reg_V_groove" || geometry == "shaft_wideU_groove"){
    if(is.null(dimensions$d)==FALSE){
      Area <- pi*((dimensions$d/2)^2)
    }
  }
  if(geometry == "shaft_cyl_tube_circ_hole" || geometry == "shaft_cyl_tube_circ_hole_pressurized"){
    if(is.null(dimensions$R)==FALSE && is.null(dimensions$r)==FALSE && is.null(dimensions$d)==FALSE){
      Area <- pi*(((dimensions$R)^2) - ((dimensions$r)^2)) - (dimensions$R - dimensions$r)*dimensions$d
    }
  }
  if(is.null(loadconditions$Sa) == FALSE){
    Snet <- loadconditions$Sa/Area
  }

  if(is.null(options$Kf)==TRUE){
    # Compute the stress concentration factor, Kt
    Kt <- stress.concentration.factor(dimensions,geometry)

    if(Kt*Snet < loadconditions$Sy){
      stop('The specimen is not expected to propagate since K_t*S_net is less than S_y')
    }

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
    }
  } else {
    Kf <- options$Kf
  }

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

  Sm <- 0.5*(stressmax + stressmin)
  stress_amp <- 0.5*(stressmax - stressmin)
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
    R = stressmin/stressmax

    straintrace <- function(reversal) strain_amp - (sig_f/E)*((0.5*(1 - R))^(1 - gam))*((reversal)^b) - eps_f*((0.5*(1 - R))^((c/b)*(1 - gam)))*((reversal)^c)
  }

  rev_2N_trace <- findzeros(straintrace, 10, 10^9)
  Life_trace <- rev_2N_trace/2

  return(list(Kf = Kf, Snet=Snet, stress_strain = stress_strain, Dstress_Dstrain = Dstress_Dstrain, stressmax = stressmax, strainmax = strainmax, stressmin = stressmin, strainmin = strainmin, strainreversals = rev_2N_trace,strainlife = Life_trace))
}
