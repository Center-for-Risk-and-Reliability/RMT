# Creep Analyzer
# Developed by Dr. Reuel Smith, 2022-2023

creep.analysis <- function(data,model=1,creepproperties = list(C = 20), units=1) {
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

  # Establish the axes labels
  Xlab <- "Time (hours)"
  Ylab <- paste(c("Stress (",stresslabel,")"),collapse = "")

  # Set up details for the parameter evaluation
  # Isolate temperatures
  temp_v <- data[,1][!duplicated(data[,1])]
  # Full temperature, stress,, and failure time vectors from data
  temp_v2 <- data[,1]
  stress_v <- data[,2]
  time_v <- data[,3]
  # Set up temperature label for plotting
  templabel <- rep(0,length(data[,1]))
  # Set time limits based on time_v
  xlims <- c(10^floor(log10(min(time_v))), 10^ceil(log10(max(time_v))))

  # Computes the upper and lower bound for the TTF axis in terms of log-time
  signs1 <- c(floor(log10(min(time_v))):ceiling(log10(max(time_v))))
  logtimes1 <- 10^signs1
  Pticks1X <- c(1:(9*length(logtimes1)-8))
  Pticks1X[1] <- logtimes1[1]
  Pticks1Xlabel <- Pticks1X

  for(i2 in 1:(length(signs1)-1)){
    Pticks1X[(9*i2-7):(9*(i2+1)-8)] <- logtimes1[i2]*c(2:10)
    Pticks1Xlabel[(9*i2-7):(9*(i2+1)-8)] <- c("","","","","","","","",logtimes1[i2+1])
  }

  for(i in 1:length(data[,1])){
    templabel[i] <- paste(c(data[,1][i]," \U00B0","F"),collapse = "")
  }

  # set up Larson-Miller, Mason-Haferd Parameter, or Sherby-Dorn P
  timelines  <- rep(0,2*length(temp_v))
  stresslines  <- rep(0,2*length(temp_v))
  bestfitlines <- rep(0,2*length(temp_v))
  P <- rep(0,length(temp_v))
  A  <- rep(0,length(temp_v))
  b  <- rep(0,length(temp_v))

  # Larson-Miller computation
  if(model == 1){
    Strace <- creepproperties$Strace
    # Pull Larson-Miller constant
    if(length(creepproperties$C) == 0){
      C = 20
    } else{
      C <- creepproperties$C
    }

    for(i in 1:length(temp_v)){
      # Set time and stress vectors for fitting
      time_fit <- time_v[which(temp_v2 == temp_v[i])]
      stress_fit <- stress_v[which(temp_v2 == temp_v[i])]

      # Compute A and b from semilog-linear fit S = A + b*log(t)
      creepparam  <- lm(stress_fit ~ poly(log(time_fit), 1, raw=TRUE))
      A[i] <- summary(creepparam)$coefficients[1,1]
      b[i] <- summary(creepparam)$coefficients[2,1]

      # Generate stress, time, and best-fit lines
      stressline_i <- A[i] + b[i]*log(c(min(time_fit),max(time_fit)))
      xlims_i <- c(min(time_fit),max(time_fit))
      stresslines[(i*2) - 1] <- stressline_i[1]
      stresslines[i*2] <- stressline_i[2]
      timelines[(i*2) - 1] <- xlims_i[1]
      timelines[i*2] <- xlims_i[2]
      bestfitlines[(i*2) - 1] <- paste(c("Best fit for ",temp_v[i]," \U00B0","F"),collapse = "")
      bestfitlines[i*2] <- paste(c("Best fit for ",temp_v[i]," \U00B0","F"),collapse = "")

      # Compute Larson-Miller coefficients
      if(length(which(stress_fit==Strace)) == 0){
        P[i] <- (temp_v[i] + Tempconv)*(C + log10(exp((Strace - A[i])/b[i])))
      }
      if(length(which(stress_fit==Strace)) == 1){
        P[i] <- (temp_v[i] + Tempconv)*(C + log10(time_fit[which(stress_fit == Strace)]))
      }
    }
    # Isolate reference temperature if any exist to label primary Larson-Miller parameter to compute equivalent test time
    if(length(creepproperties$Tref)==1){
      Pref <- P[which(temp_v == creepproperties$Tref)]
      time_test <- 10^(Pref/(creepproperties$Tref + Tempconv) - C)
      time_test_v <- 10^(Pref/(temp_v + Tempconv) - C)
    }
  }

  # Mason-Haferd computation
  if(model == 2){
    Strace <- creepproperties$Strace
    # Pull Mason-Haferd constants
    temp_a <- creepproperties$temp_a
    log10t_a <- creepproperties$log10t_a

    for(i in 1:length(temp_v)){
      # Set time and stress vectors for fitting
      time_fit <- time_v[which(temp_v2 == temp_v[i])]
      stress_fit <- stress_v[which(temp_v2 == temp_v[i])]

      # Compute A and b from semilog-linear fit S = A + b*log(t)
      creepparam  <- lm(stress_fit ~ poly(log(time_fit), 1, raw=TRUE))
      A[i] <- summary(creepparam)$coefficients[1,1]
      b[i] <- summary(creepparam)$coefficients[2,1]

      # Generate stress, time, and best-fit lines
      stressline_i <- A[i] + b[i]*log(c(min(time_fit),max(time_fit)))
      xlims_i <- c(min(time_fit),max(time_fit))
      stresslines[(i*2) - 1] <- stressline_i[1]
      stresslines[i*2] <- stressline_i[2]
      timelines[(i*2) - 1] <- xlims_i[1]
      timelines[i*2] <- xlims_i[2]
      bestfitlines[(i*2) - 1] <- paste(c("Best fit for ",temp_v[i]," \U00B0","F"),collapse = "")
      bestfitlines[i*2] <- paste(c("Best fit for ",temp_v[i]," \U00B0","F"),collapse = "")

      # Compute Mason-Haferd coefficients
      if(length(which(stress_fit==Strace)) == 0){
        P[i] <- (temp_v[i] - temp_a)/(log10(exp((Strace - A[i])/b[i])) - log10t_a)
      }
      if(length(which(stress_fit==Strace)) == 1){
        P[i] <- (temp_v[i] - temp_a)/(log10(time_fit[which(stress_fit == Strace)])  - log10t_a)
      }
    }
    # Isolate reference temperature if any exist to label primary Mason-Haferd parameter to compute equivalent test time
    if(length(creepproperties$Tref)==1){
      Pref <- P[which(temp_v == creepproperties$Tref)]
      time_test <- 10^((creepproperties$Tref - temp_a)/Pref + log10t_a)
      time_test_v <- 10^((temp_v - temp_a)/Pref + log10t_a)
    }
  }

  # Sherby-Dorn computation
  if(model == 3){
    Strace <- creepproperties$Strace
    # Pull Activation Energy
    if(length(creepproperties$Q)>0){
      Q <- creepproperties$Q
    }
    if(length(creepproperties$Ea)>0){
      Ea <- creepproperties$Ea
    }

    # Convert temperature as necessary
    if(units == 1){
      temp_v_SherbyDorn <- temp_v + 273.15
    }
    if(units == 2){
      temp_v_SherbyDorn <- (temp_v + 459.67)*(5/9)
    }
    R <- 8.314
    k_B <- 8.617e-5

    for(i in 1:length(temp_v)){
      # Set time and stress vectors for fitting
      time_fit <- time_v[which(temp_v2 == temp_v[i])]
      stress_fit <- stress_v[which(temp_v2 == temp_v[i])]

      # Compute A and b from semilog-linear fit S = A + b*log(t)
      creepparam  <- lm(stress_fit ~ poly(log(time_fit), 1, raw=TRUE))
      A[i] <- summary(creepparam)$coefficients[1,1]
      b[i] <- summary(creepparam)$coefficients[2,1]

      # Generate stress, time, and best-fit lines
      stressline_i <- A[i] + b[i]*log(c(min(time_fit),max(time_fit)))
      xlims_i <- c(min(time_fit),max(time_fit))
      stresslines[(i*2) - 1] <- stressline_i[1]
      stresslines[i*2] <- stressline_i[2]
      timelines[(i*2) - 1] <- xlims_i[1]
      timelines[i*2] <- xlims_i[2]
      bestfitlines[(i*2) - 1] <- paste(c("Best fit for ",temp_v[i]," \U00B0","F"),collapse = "")
      bestfitlines[i*2] <- paste(c("Best fit for ",temp_v[i]," \U00B0","F"),collapse = "")

      # Compute Sherby-Dorn coefficients
      if(length(which(stress_fit==Strace)) == 0){
        if(length(creepproperties$Q)>0){
          P[i] <- log10(exp((Strace - A[i])/b[i])) - 0.43*(Q/(R*temp_v_SherbyDorn[i]))
        }
        if(length(creepproperties$Ea)>0){
          P[i] <- log10(exp((Strace - A[i])/b[i])) - 0.43*(Ea/(k_B*temp_v_SherbyDorn[i]))
        }
      }
      if(length(which(stress_fit==Strace)) == 1){
        if(length(creepproperties$Q)>0){
          P[i] <- log10(time_fit[which(stress_fit == Strace)]) - 0.43*(Q/(R*temp_v_SherbyDorn[i]))
        }
        if(length(creepproperties$Ea)>0){
          P[i] <- log10(time_fit[which(stress_fit == Strace)]) - 0.43*(Ea/(k_B*temp_v_SherbyDorn[i]))
        }
      }
    }
    # Isolate reference temperature if any exist to label primary Larson-Miller parameter to compute equivalent test time
    if(length(creepproperties$Tref)==1){
      if(units == 1){
        Tref_SherbyDorn <- creepproperties$Tref + 273.15
      }
      if(units == 2){
        Tref_SherbyDorn <- (creepproperties$Tref + 459.67)*(5/9)
      }
      Pref <- P[which(temp_v == creepproperties$Tref)]
      if(length(creepproperties$Q)>0){
        time_test <- 10^(Pref + 0.43*(Q/(R*Tref_SherbyDorn)))
        time_test_v <- 10^(Pref + 0.43*(Q/(R*temp_v_SherbyDorn)))
      }
      if(length(creepproperties$Ea)>0){
        time_test <- 10^(Pref + 0.43*(Ea/(k_B*Tref_SherbyDorn)))
        time_test_v <- 10^(Pref + 0.43*(Ea/(k_B*temp_v_SherbyDorn)))
      }
    }
  }
  df <- data.frame(time = data[,3], stress = data[,2], Temperature = templabel)
  df2 <- data.frame(time2 = timelines, stress2 = stresslines, Best_Fit = bestfitlines)

  plotout<-ggplot() +
    geom_point(data=df, aes(time,stress, shape = Temperature), colour = 'black', size = 2.2) +
    scale_x_continuous(trans = 'log10', limits = xlims, breaks=Pticks1X, labels=Pticks1Xlabel) +
    scale_y_continuous(limits = c(0,max(stresslines))) +
    labs(x=Xlab,y=Ylab)
  plotout <- plotout + geom_path(data=df2, aes(time2,stress2, colour = Best_Fit), size = 0.9, linetype = "solid")

  if(length(creepproperties$Tref)==1){
    return(list(CreepPlot=plotout,TestTemp=temp_v,RefCreepModelParam=Pref,RuptureTimeatRefCreepModelParam=time_test_v))
  } else{
    return(list(CreepPlot=plotout,TestTemp=temp_v,CreepModelParam=P))
  }

}
