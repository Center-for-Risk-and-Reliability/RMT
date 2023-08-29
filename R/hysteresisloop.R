# Hysteresis Loop Plot
# Developed by Reuel Smith, 2022

hysteresisloop.plot <- function(E,K,n,stressunits,loadconditions){
  library(pracma)
  library(nls.multstart)
  library(ggplot2)

  if(missing(E)){
    stop('Enter modulus of elasticity E (in MPa or ksi).')
  }

  # Check units and set up axis labels
  if(missing(stressunits) || stressunits == 1){
    stressunits <- c("MPa")
  }
  if(stressunits == 2){
    stressunits <- c("ksi")
  }

  # Check units and set up axis labels
  Xlab <- paste(c("Total Strain, (eps_tot)"),collapse="", sep = "_")
  Ylab <- paste(c("Alternating Stress, S_a (",stressunits,")"),collapse="", sep = "_")

  # Load conditions (loadconditions) may be made up of single min or max stress or strain.
  # It is assumed to be fully reversed if no other bound is given.
  if(length(loadconditions$maxstress) == 1 && length(loadconditions$minstress) == 0){
    maxstress <- loadconditions$maxstress
    minstress <- -maxstress
    Dstress <- 2*maxstress
  }
  if(length(loadconditions$maxstress) == 0 && length(loadconditions$minstress) == 1){
    minstress <- loadconditions$minstress
    maxstress <- -minstress
    Dstress <- 2*maxstress
  }
  if(length(loadconditions$maxstrain) == 1 && length(loadconditions$minstrain) == 0){
    maxstrain <- loadconditions$maxstrain
    minstrain <- -maxstrain
    Dstrain <- 2*maxstrain
  }
  if(length(loadconditions$maxstrain) == 0 && length(loadconditions$minstrain) == 1){
    minstrain <- loadconditions$minstrain
    maxstrain <- -minstrain
    Dstrain <- 2*maxstrain
  }
  # Load conditions (loadconditions) may also be made up of a stress range or a strain range
  if(length(loadconditions$stressrange) == 2 && length(loadconditions$strainrange) == 0){
    maxstress <- max(loadconditions$stressrange)
    minstress <- min(loadconditions$stressrange)
    Dstress <- maxstress - minstress
  }
  if(length(loadconditions$stressrange) == 0 && length(loadconditions$strainrange) == 2){
    maxstrain <- max(loadconditions$strainrange)
    minstrain <- min(loadconditions$strainrange)
    Dstrain <- maxstrain - minstrain
  }

  # Start with the initial loading.  Assume tension (loadtype = 1) first unless stated.
  if(length(loadconditions$loadtype) == 0 || (length(loadconditions$loadtype) == 1 && loadconditions$loadtype == 1)){
    loadtype <- 1
  }
  # Otherwise it is compression (loadtype = 2)
  if(length(loadconditions$loadtype) == 1 && loadconditions$loadtype == 2){
    loadtype <- 2
  }
  # Now set up initial loading curve based on stress
  if(length(loadconditions$maxstress) == 1 || length(loadconditions$minstress) == 1 || length(loadconditions$stressrange) == 2){
    if(loadtype == 1){
      stressline1 <- linspace(0,maxstress,100)
      stressline2 <- linspace(maxstress,minstress,100)
      Dstress2 <- maxstress - stressline2
      stressline3 <- linspace(minstress,maxstress,100)
      Dstress3 <- stressline3 - minstress
    }
    if(loadtype == 2){
      stressline1 <- linspace(0,minstress,100)
      stressline2 <- linspace(minstress,maxstress,100)
      Dstress2 <- stressline2 - minstress
      stressline3 <- linspace(maxstress,minstress,100)
      Dstress3 <- maxstress - stressline3
    }
    strainline1 <- (stressline1/E) + (stressline1/K)^(1/n)
    Dstrain2 <- (Dstress2/E) + 2*((Dstress2/(2*K))^(1/n))
    Dstrain3 <- (Dstress3/E) + 2*((Dstress3/(2*K))^(1/n))
    if(loadtype == 1){
      strainline2 <- max(strainline1) - Dstrain2
      strainline3 <- min(strainline2) + Dstrain3
    }
    if(loadtype == 2){
      strainline2 <- min(strainline1) + Dstrain2
      strainline3 <- max(strainline2) - Dstrain3
    }
  }

  # Now set up initial loading curve based on strain
  if(length(loadconditions$maxstrain) == 1 || length(loadconditions$minstrain) == 1 || length(loadconditions$strainrange) == 2){
    if(loadtype == 1){
      strainline1 <- linspace(0,maxstrain,100)
      strainline2 <- linspace(maxstrain,minstrain,100)
      Dstrain2 <- maxstrain - strainline2
      strainline3 <- linspace(minstrain,maxstrain,100)
      Dstrain3 <- strainline3 - minstrain
    }
    stressline1 <- rep(0,100)
    Dstress2 <- rep(0,100)
    Dstress3 <- rep(0,100)

    for(i in 1:100){
      hysloop1 <- function(stress) strainline1[i] - (stress/E) - ((stress/K)^(1/n))
      stressline1[i] <- findzeros(hysloop1, 0, K)
      hysloop2 <- function(Dstress) Dstrain2[i] - (Dstress/E) - 2*((Dstress/(2*K))^(1/n))
      Dstress2[i] <- findzeros(hysloop2, 0, 2*K)
      hysloop3 <- function(Dstress) Dstrain3[i] - (Dstress/E) - 2*((Dstress/(2*K))^(1/n))
      Dstress3[i] <- findzeros(hysloop3, 0, 2*K)
    }
    stressline2 <- max(stressline1) - Dstress2
    stressline3 <- min(stressline2) + Dstress3
  }

  maxstrain <- max(strainline2)
  minstrain <- min(strainline2)
  maxstress <- max(stressline2)
  minstress <- min(stressline2)

  # =======================================
  # Plotting Output
  # =======================================
  df <- data.frame(strain = c(strainline1,strainline2,strainline3), stress = c(stressline1,stressline2,stressline3), curvelines = c(rep("Initial Loading",100),rep("Hysteresis Loop from Repeated Loading",200)))

  plotout1<-ggplot() +
    geom_path(data=df, aes(strain, stress, colour = curvelines), size = 0.9) +
    xlab(Xlab) +
    ylab(Ylab)

  return(list(hysteresisloop = plotout1,maxstrain = maxstrain,minstrain = minstrain,maxstress = maxstress,minstress = minstress))
}
