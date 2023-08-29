# Variable Amplitude Loading Modeling with Damage
# Developed by Reuel Smith, 2022

var.amp.loadingdamage.model <- function(dat,damagerule="Miner",stressunits){
  # dat is entered as a list made up of Nf, ni, Di, stress blocks,
  # corr_rel (correction relationship), reversals at stress blocks, sig_f, and/or b
  library(pracma)
  library(ggplot2)

  if(missing(dat)){
    stop('Enter data (in list form).')
  }
  if(is.list(dat)==FALSE){
    stop('Enter data as list.')
  }

  # Check units and set up axis labels
  if(missing(stressunits) || stressunits == 1){
    stressunits <- c("MPa")
  }
  if(stressunits == 2){
    stressunits <- c("ksi")
  }

  # Check damage rule (damagerule) conditions (1 for linear: Palmgren-Miner,
  # 2 for non-linear: Kwofie & Rahbar, 3 for non-linear: Corten-Dolan)
  # Add conditions for non-linear at a later date and add it back into the input as well

  # Computation 1: Operation cycles (ni) from cycles to failure (Nfi) and damage percentage (Dfi)
  if(length(dat$Nf) > 1 && length(dat$Df) > 1 && length(dat$Nf) != length(dat$Df)){
    stop('Check to see that your cycle to failure vector (Nf) and damage vector (Df) are the same length')
  }
  if(length(dat$Nf) > 1 && length(dat$Df) > 1 && length(dat$Nf) == length(dat$Df)){
    # Check damage sum
    if(sum(dat$Df) != 1){
      stop('Make sure the sum of your damage is equal to 1.')
    }
    if(damagerule == "Miner"){
      # Calculate total operation cycles
      n <- 1/sum(dat$Df/dat$Nf)
      ni <- n*dat$Df
    }
    if(damagerule == "KwofieRahbar"){
      # Calculate total operation cycles
      n <- 1/sum((dat$Df/dat$Nf)*(log(dat$Nf)/log(dat$Nf[1])))
      ni <- n*dat$Df
    }
    return(list(opercyclesbyblock = ni, totalopercycles = n))
  }

  # Computation 2: Blocks and damage percentage from operation cycles (ni) and cycles to failure (Nfi).
  if(length(dat$Nf) > 1 && length(dat$ni) > 1 && length(dat$Nf) != length(dat$ni)){
    stop('Check to see that your cycle to failure vector (Nf) and operation cycles vector (ni) are the same length')
  }
  if(length(dat$Nf) > 1 && length(dat$ni) > 1 && length(dat$Nf) == length(dat$ni)){
    if(damagerule == "Miner"){
      # Calculate total operation damage and blocks to failure
      Bfi <- dat$ni/dat$Nf
      Bf <- 1/sum(dat$ni/dat$Nf)
      Di <- Bf*(dat$ni/dat$Nf)
    }
    if(damagerule == "KwofieRahbar"){
      # Calculate total operation damage and blocks to failure
      Bfi <- (dat$ni/dat$Nf)*(log(dat$Nf)/log(Nf[1]))
      Bf <- 1/sum((dat$ni/dat$Nf)*(log(dat$Nf)/log(dat$Nf[1])))
      Di <- Bf*((dat$ni/dat$Nf)*(log(dat$Nf)/log(dat$Nf[1])))
    }
    return(list(damagebyblock = Di, repetitionstofailure = Bf))
  }

  # Computation 3: Operation cycles (ni) from sig_f, b, stress ranges (sranges), and damage percentage (Dfi)
  if(length(dat$sranges) > 1 && length(dat$Df) > 1 && length(dat$sig_f) == 1 && length(dat$b) == 1 && length(dat$sranges) != length(dat$Df)){
    stop('Check to see that your stress range list (sranges) and damage vector (Df) are the same length')
  }
  if(length(dat$sranges) > 1 && length(dat$Df) > 1 && length(dat$sig_f) == 1 && length(dat$b) == 1 && length(dat$sranges) == length(dat$Df)){
    # Check damage sum
    if(sum(dat$Df) != 1){
      stop('Make sure the sum of your damage is equal to 1.')
    }
    # Now compute Nfs (may be separated soon)
    Nf<-rep(0,length(dat$sranges))
    Sar_m<-rep(0,length(dat$sranges))
    for(i in 1:length(dat$Df)){
      # Stress amplitude
      Sa <- 0.5*(max(dat$sranges[[i]]) - min(dat$sranges[[i]]))
      Sm <- 0.5*(max(dat$sranges[[i]]) + min(dat$sranges[[i]]))
      if(Sm == 0){
        Sar <- Sa
      } else{
        # Pull Mean Stress Correction
        # Check correction relationships dat$corr_rel and needed input
        if(length(is.na(dat$corr_rel))==0){
          stop('Please enter a correction relationship as "corr_rel".')
        } else{
          if(dat$corr_rel == "Soderberg" && length(dat$Sy) == 1){
            Sy <- dat$Sy
            Sar <- Sa/(1 - (Sm/dat$Sy))
          }
          if(dat$corr_rel == "ModifiedGoodman" && length(Su) == 1){
            Sar <- Sa/(1 - (Sm/dat$Su))
          }
          if(dat$corr_rel == "Morrow"){
            Sar <- Sa/(1 - (Sm/dat$sig_f))
          }
          if(dat$corr_rel == "Gerber" && length(Su) == 1){
            Sar <- Sa/(1 - ((Sm/dat$Su)^2))
          }
          if(dat$corr_rel == "SWT"){
            # Stress Ratio
            if(max(dat$sranges[[i]]) == 0){
              stop('Select another correlation relationship model other than SWT or Walker.')
            }
            R <- min(dat$sranges[[i]])/max(dat$sranges[[i]])
            Sar <- Sa/sqrt(0.5*(1 - R))
          }
          if(dat$corr_rel == "Walker" && length(dat$gam) == 1){
            # Stress Ratio
            if(max(dat$sranges[[i]]) == 0){
              stop('Select another correlation relationship model other than SWT or Walker.')
            }
            R <- min(dat$sranges[[i]])/max(dat$sranges[[i]])
            Sar <- Sa/((0.5*(1 - R))^(1-dat$gam))
          }
        }
      }
      stresstrace <- function(reversal) Sar - dat$sig_f*(reversal)^dat$b
      rev <- findzeros(stresstrace, 10, 10^9)
      Nf[i] <- rev/2
    }
    if(damagerule == "Miner"){
      # Calculate total operation cycles
      n <- 1/sum(dat$Df/Nf)
      ni <- n*dat$Df
    }
    if(damagerule == "KwofieRahbar"){
      # Calculate total operation cycles
      n <- 1/sum((dat$Df/Nf)*(log(dat$Nf)/log(Nf[1])))
      ni <- n*dat$Df
    }
    return(list(opercyclesbyblock = ni, cyclestofailurebyblock = Nf, reversalstofailurebyblock = 2*Nf, totalopercycles = n))
  }
  # Computation 4: Blocks and damage percentage from sig_f, b, stress ranges (sranges), and operation cycles (ni).
  if(length(dat$sranges) > 1 && length(dat$ni) > 1 && length(dat$sig_f) == 1 && length(dat$b) == 1 && length(dat$sranges) != length(dat$ni)){
    stop('Check to see that your cycle to failure vector (Nf) and operation cycles vector (ni) are the same length')
  }
  # Now compute Nfs (may be separated soon)
  Nf<-rep(0,length(dat$sranges))
  Sar_m<-rep(0,length(dat$sranges))
  for(i in 1:length(dat$ni)){
    # Stress amplitude
    Sa <- 0.5*(max(dat$sranges[[i]]) - min(dat$sranges[[i]]))
    Sm <- 0.5*(max(dat$sranges[[i]]) + min(dat$sranges[[i]]))
    if(Sm == 0){
      Sar <- Sa
    } else{
      # Pull Mean Stress Correction
      # Check correction relationships dat$corr_rel and needed input
      if(length(is.na(dat$corr_rel))==0){
        stop('Please enter a correction relationship as "corr_rel".')
      } else{
        if(dat$corr_rel == "Soderberg" && length(dat$Sy) == 1){
          Sy <- dat$Sy
          Sar <- Sa/(1 - (Sm/dat$Sy))
        }
        if(dat$corr_rel == "ModifiedGoodman" && length(Su) == 1){
          Sar <- Sa/(1 - (Sm/dat$Su))
        }
        if(dat$corr_rel == "Morrow"){
          Sar <- Sa/(1 - (Sm/dat$sig_f))
        }
        if(dat$corr_rel == "Gerber" && length(Su) == 1){
          Sar <- Sa/(1 - ((Sm/dat$Su)^2))
        }
        if(dat$corr_rel == "SWT"){
          # Stress Ratio
          if(max(dat$sranges[[i]]) == 0){
            stop('Select another correlation relationship model other than SWT or Walker.')
          }
          R <- min(dat$sranges[[i]])/max(dat$sranges[[i]])
          Sar <- Sa/sqrt(0.5*(1 - R))
        }
        if(dat$corr_rel == "Walker" && length(dat$gam) == 1){
          # Stress Ratio
          if(max(dat$sranges[[i]]) == 0){
            stop('Select another correlation relationship model other than SWT or Walker.')
          }
          R <- min(dat$sranges[[i]])/max(dat$sranges[[i]])
          Sar <- Sa/((0.5*(1 - R))^(1-dat$gam))
        }
      }
    }
    stresstrace <- function(reversal) Sar - dat$sig_f*((reversal)^dat$b)
    rev <- findzeros(stresstrace, 10, 10^9)
    Nf[i] <- rev/2
    Sar_m[i]<-Sar
  }
  if(length(dat$sranges) > 1 && length(dat$ni) > 1 && length(dat$sig_f) == 1 && length(dat$b) == 1 && length(dat$sranges) == length(dat$ni)){
    if(damagerule == "Miner"){
      # Calculate total operation damage and blocks to failure
      Bfi <- dat$ni/Nf
      Bf <- 1/sum(dat$ni/Nf)
      Di <- Bf*(dat$ni/Nf)
    }
    if(damagerule == "KwofieRahbar"){
      # Calculate total operation damage and blocks to failure
      Bfi <- (dat$ni/Nf)*(log(Nf)/log(Nf[1]))
      Bf <- 1/sum((dat$ni/Nf)*(log(Nf)/log(Nf[1])))
      Di <- Bf*(dat$ni/Nf)*(log(Nf)/log(Nf[1]))
    }
    return(list(damagebyblock = Di, cyclestofailurebyblock = Nf, reversalstofailurebyblock = 2*Nf, repetitionstofailure = Bf))
  }
}
