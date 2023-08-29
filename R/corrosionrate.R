# Corrosion Rate Calculator
# Developed by Dr. Reuel Smith, 2022

corrosion.rate <- function(data,model=1,creepproperties = list(C = 20), units=1) {
  library(ggplot2)
  library(pracma)

  # Check units and set up axis labels
  if(units == 1){
    stresslabel <- c("MPa")
    Tempconv <- 273.15
  }
  if(units == 2){
    stresslabel <- c("ksi")
    Tempconv <- 459.67
  }

  # Stern-Gary Model: i_corr = B/R_p (log i_corr = log B - log R_p) A/cm^2
  # Yalcyn-Ergun Model: i_corr = i_0 exp(-Ct) (log i_corr = log i_0 - Ct)
  # Liu-Weyers' Model: i_corr = 0.92 exp(8.37 + 0.618 log(1.69 CL) - 3034/T - 0.000105 R_c +2.32t^-0.215)
  # (log i_corr = log 0.92 + 8.37 + 0.618 log(1.69 CL) - 3034/T - 0.000105 R_c +2.32t^-0.215) muA/cm^2
  # Thickness Loss: i_corr= (j x m)/(n x F x rho) cm/s
  #                 i_corr= 0.00327 x (j x m)/(n x rho) mm/year
}
