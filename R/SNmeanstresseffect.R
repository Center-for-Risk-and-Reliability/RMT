# Stress-Life Mean Stress Effect Calculator
# Developed by Reuel Smith, 2022

SN.meanstresseffect <- function(relationship,Sa,Sm,altvar){
  if(relationship == "Soderberg"){
    Sy <- altvar$Sy
    Sar <- Sa/(1 - (Sm/Sy))
  }
  if(relationship == "ModGoodman"){
    Su <- altvar$Su
    Sar <- Sa/(1 - (Sm/Su))
  }
  if(relationship == "Morrow"){
    sig_f <- altvar$sig_f
    Sar <- Sa/(1 - (Sm/sig_f))
  }
  if(relationship == "Gerber"){
    Su <- altvar$Su
    Sar <- Sa/(1 - ((Sm/Su)^2))
  }
  if(relationship == "SWT"){
    R <- altvar$R
    Sar <- Sa/sqrt(0.5*(1-R))
  }
  if(relationship == "Walker"){
    R <- altvar$R
    gam <- altvar$gam
    Sar <- Sa/((0.5*(1-R))^(1-gam))
  }
  return(Sar)
}
