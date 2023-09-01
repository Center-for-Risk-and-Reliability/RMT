# Pitting Corrosion Calculator
# Developed by Dr. Reuel Smith, 2022-2023

corrosion.pitting <- function(data, r_0 = 0, r_cr = 1, matproperties = NULL, units=1) {
  library(ggplot2)
  library(pracma)
  library(stringr)

  # Check units and set up axis labels
  if(units == 1){
    colnames(data)[4] <- "Temperature (K)"
    Tempconv <- 273.15
  }
  if(units == 2){
    colnames(data)[4] <- "Temperature (K)"
    data[,4] <- data[,4]*(5/9)
  }
  # Boltzmann Constant (eV/K)
  K <- 8.617385e-5

  # Kondo-Wei Model: i_corr = 2 pi r^2 dr/dt = ((M I_p0)/(n F rho)) exp(-E_a/kT)
  #                  Dt = ((2 pi)/3)((n F rho)/(M I_p0))(r_f^3 - r_0^3) exp(E_a/kT)
  # OR
  #                  i_corr = ((2 pi)/(3 A pH)) exp(-E_a/kT)
  #                  Dt = A pH (r_f^3 - r_0^3) exp(E_a/kT)
  # First model would seek I_p0 and E_a as parameters while the second would seek A and E_a
  # Data may be used directly with ADT tools once Kondo-Wei model is implemented.

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # Pulls names from column headers for the data
  col_names_data <- colnames(data)

  time_output_names <- c(col_names_data[1],"Censored",col_names_data[4:length(col_names_data)])

  # Pulls stress from columns 4 and up
  stresscount<-dim(data)[2]-3
  if(stresscount>1){
    stressvals<-vector(mode = "list", length = stresscount)
  }

  for(i in 1:length(unitnames)){
    stressgroup<-which(data[,3]==unitnames[i])
    stressvals0<-data[stressgroup[1],4:dim(data)[2]]
    names(stressvals0)<-NULL
    if(stresscount==1){
      if(i==1){
        stressvals<-stressvals0
      } else if(i>1){
        stressvals<-c(stressvals,stressvals0)
      }
    } else {
      for(j in 1:stresscount){
        if(i==1){
          stressvals[[j]]<-stressvals0[j]
        } else if(i>1){
          stressvals[[j]]<-c(stressvals[[j]],stressvals0[j])
        }
      }
    }
  }

  # Preliminary transformation of units in the event of material properties that would necessitate the full
  # expression of Kondo-Wei
  if(length(matproperties) > 0){
    # Check time for if it is in hours or minutes (it needs to be in seconds)
    if(isTRUE(str_detect(col_names_data[1], "hr", negate = FALSE)) || isTRUE(str_detect(col_names_data[1], "hour", negate = FALSE))){
      data[,1] <- 3600*data[,1]
    }
    if(isTRUE(str_detect(col_names_data[1], "min", negate = FALSE))){
      data[,1] <- 60*data[,1]
    }
    # Check radius if it is in mm or inches (it needs to be in meters)
    if(isTRUE(str_detect(col_names_data[2], "mm", negate = FALSE))){
      data[,2] <- (1/1000)*data[,2]
      r_0 <- (1/1000)*r_0
      r_cr <- (1/1000)*r_cr
    }
    if(isTRUE(str_detect(col_names_data[2], "in", negate = FALSE)) || isTRUE(str_detect(col_names_data[2], "inch", negate = FALSE))){
      data[,2] <- 0.0254*data[,2]
      r_0 <- 0.0254*r_0
      r_cr <- 0.0254*r_cr
    }
    time_output_names[1] <- "Time (seconds)"
    xlabel1 <- "Time (seconds)"
    ylabel1 <- "Pit Radius (meters)"
    Alabel1 <- "A (sec/m\U00B3)"
    icorrlabel1 <- "i_corr (m\U00B3/sec)"
  }
  if(length(matproperties) == 0){
    if(isTRUE(str_detect(col_names_data[1], "hr", negate = FALSE)) || isTRUE(str_detect(col_names_data[1], "hour", negate = FALSE))){
      time_output_names[1] <- "Time (hours)"
      xlabel1 <- "Time (hours)"
    }
    if(isTRUE(str_detect(col_names_data[1], "min", negate = FALSE))){
      time_output_names[1] <- "Time (minutes)"
      xlabel1 <- "Time (minutes)"
    }
    # Check radius if it is in mm, meters, or inches
    if(isTRUE(str_detect(col_names_data[2], "mm", negate = FALSE))){
      ylabel1 <- "Pit Radius (mm)"
    }
    if(isTRUE(str_detect(col_names_data[2], "m", negate = FALSE)) && isFALSE(str_detect(col_names_data[2], "mm", negate = FALSE))){
      ylabel1 <- "Pit Radius (m)"
    }
    if(isTRUE(str_detect(col_names_data[2], "in", negate = FALSE)) || isTRUE(str_detect(col_names_data[2], "inch", negate = FALSE))){
      ylabel1 <- "Pit Radius (inches)"
    }
    if((isTRUE(str_detect(col_names_data[1], "hr", negate = FALSE)) || isTRUE(str_detect(col_names_data[1], "hour", negate = FALSE))) && isTRUE(str_detect(col_names_data[2], "mm", negate = FALSE))){
      Alabel1 <- "A (hour/mm\U00B3)"
      icorrlabel1 <- "i_corr (mm\U00B3/hour)"
    }
    if(isTRUE(str_detect(col_names_data[1], "min", negate = FALSE)) && isTRUE(str_detect(col_names_data[2], "mm", negate = FALSE))){
      Alabel1 <- "A (min/mm\U00B3)"
      icorrlabel1 <- "i_corr (mm\U00B3/min)"
    }
    if((isTRUE(str_detect(col_names_data[1], "hr", negate = FALSE)) || isTRUE(str_detect(col_names_data[1], "hour", negate = FALSE))) && (isTRUE(str_detect(col_names_data[2], "m", negate = FALSE)) && isFALSE(str_detect(col_names_data[2], "mm", negate = FALSE)))){
      Alabel1 <- "A (hour/m\U00B3)"
      icorrlabel1 <- "i_corr (m\U00B3/hour)"
    }
    if(isTRUE(str_detect(col_names_data[1], "min", negate = FALSE)) && (isTRUE(str_detect(col_names_data[2], "m", negate = FALSE)) && isFALSE(str_detect(col_names_data[2], "mm", negate = FALSE)))){
      Alabel1 <- "A (min/m\U00B3)"
      icorrlabel1 <- "i_corr (m\U00B3/min)"
    }
    if((isTRUE(str_detect(col_names_data[1], "hr", negate = FALSE)) || isTRUE(str_detect(col_names_data[1], "hour", negate = FALSE))) && (isTRUE(str_detect(col_names_data[2], "in", negate = FALSE)) || isTRUE(str_detect(col_names_data[2], "inch", negate = FALSE)))){
      Alabel1 <- "A (hour/in\U00B3)"
      icorrlabel1 <- "i_corr (in\U00B3/hour)"
    }
    if(isTRUE(str_detect(col_names_data[1], "min", negate = FALSE)) && (isTRUE(str_detect(col_names_data[2], "in", negate = FALSE)) || isTRUE(str_detect(col_names_data[2], "inch", negate = FALSE)))){
      Alabel1 <- "A (min/in\U00B3)"
      icorrlabel1 <- "i_corr (in\U00B3/min)"
    }
  }

  # Case 1: Input data table (t, r, T) and material properties M (g/mol), n, rho (g/m^3)
  #         Output t vs. r plot, t vs. volume plot, I_p0 and E_a
  if (dim(data)[2] == 4 && length(matproperties) > 0){
    # r^3 = r_0^3 + (Dt/(((2 pi)/3)((n F rho)/(M I_p0)))) exp(-E_a/kT)
    # theta[1] ~ I_p0, theta[2] ~ E_a

    # Faraday Constant in C/mol
    Farad_C <-96514

    lifedamoutput <- function(Lifedat,Damdat,Tempdat,Dam_fail){
      # r_0 and Damdat in meters
      params  <- pinv(matrix(c(rep(1,length(Damdat)),-1/Tempdat),nrow = length(Damdat), ncol = 2, byrow = FALSE))%*%(log((Damdat^3)-(r_0^3)) - log(Lifedat) + log((2/3)*pi*(matproperties$n*Farad_C*matproperties$rho)*(1/matproperties$M)))
      lifedamparams <- c(exp(params[1]),params[2]*K)
      timepsuedo <- ((2/3)*pi*(matproperties$n*Farad_C*matproperties$rho)*(1/(exp(params[1])*matproperties$M)))*((r_cr^3) - (r_0^3))*exp(params[2]*(1/Tempdat[1]))
      icorr <- ((matproperties$M*lifedamparams[1])/(matproperties$n*matproperties$rho*Farad_C))*exp(-lifedamparams[2]/(K*Tempdat[1]))

      Dest <- ((r_0^3) + ((3*lifedamparams[1]*matproperties$M*Lifedat)/(2*pi*matproperties$n*matproperties$rho*Farad_C))*exp(-lifedamparams[2]/(K*Tempdat)))^(1/3)

      R2 <- 1 - sum(((Damdat^3) - (Dest^3))^2)/sum(((Damdat^3) - mean(Damdat^3))^2)
      return(list(lifedamparams,timepsuedo,icorr,R2))
    }
    damfit <- function(TimeDamfit,Temp1,params){
      # TimeDamfit in seconds
      # Temp1 in Kelvin
      # r_0 in meters
      ((r_0^3) + ((3*params[1]*matproperties$M*TimeDamfit)/(2*pi*matproperties$n*matproperties$rho*Farad_C))*exp(-params[2]/(K*Temp1)))^(1/3)
    }
    table1labels <- c("I_p0 (C/s)","E_a (eV)",time_output_names[1],icorrlabel1,"R\U00B2")
  }

  # Case 2: Input data table (t, r, T, pH)
  #         Output t vs. r plot, t vs. volume plot, A and E_a
  if (dim(data)[2] == 5 && length(matproperties) == 0){
    # r^3 = r_0^3 + (Dt/(A x pH)) exp(-E_a/kT)
    # theta[1] ~ A, theta[2] ~ E_a
    lifedamoutput <- function(Lifedat,Damdat,Tempdat,pHdat,Dam_fail){
      params  <- pinv(matrix(c(rep(-1,length(Damdat)),-1/Tempdat),nrow = length(Damdat), ncol = 2, byrow = FALSE))%*%(log((Damdat^3)-(r_0^3)) - log(Lifedat) + log(pHdat))
      lifedamparams <- c(exp(params[1]),params[2]*K)
      timepsuedo <- exp(params[1])*((r_cr^3) - (r_0^3))*pHdat[1]*exp(params[2]*(1/Tempdat[1]))
      icorr <- ((2*pi)/(3*lifedamparams[1]*pHdat[1]))*exp(-lifedamparams[2]/(K*Tempdat[1]))

      Dest <- ((r_0^3) + (Lifedat/(lifedamparams[1]*pHdat))*exp(-lifedamparams[2]/(K*Tempdat)))^(1/3)

      R2 <- 1 - sum(((Damdat^3) - (Dest^3))^2)/sum(((Damdat^3) - mean(Damdat^3))^2)
      return(list(lifedamparams,timepsuedo,icorr,R2))
    }
    damfit <- function(TimeDamfit,pH1,Temp1,params){
      ((r_0^3) + (TimeDamfit/(params[1]*pH1))*exp(-params[2]/(K*Temp1)))^(1/3)
    }
    table1labels <- c(Alabel1,"E_a (eV)",time_output_names[1],icorrlabel1,"R\U00B2")
  }

  # Case 3: Input data table (t, r, T, pH) and material properties M,n, and rho
  #         Output t vs. r plot, t vs. volume plot, I_p0, A and E_a
  if (dim(data)[2] == 5 && length(matproperties) > 0){
    # r^3 = r_0^3 + (Dt/(A x pH)) exp(-E_a/kT)
    # theta[1] ~ A, theta[2] ~ I_p0, theta[3] ~ E_a

    # Faraday Constant in C/mol
    Farad_C <-96514

    lifedamoutput <- function(Lifedat,Damdat,Tempdat,pHdat,Dam_fail){
      params  <- pinv(matrix(c(rep(-1,length(Damdat)),-1/Tempdat),nrow = length(Damdat), ncol = 2, byrow = FALSE))%*%(log((Damdat^3)-(r_0^3)) - log(Lifedat) + log(pHdat))
      lifedamparams <- c(exp(params[1]),(2/3)*pi*((matproperties$n*Farad_C*matproperties$rho)/(matproperties$M*pHdat[1]*exp(params[1]))),params[2]*K)
      timepsuedo <- exp(params[1])*((r_cr^3) - (r_0^3))*pHdat[1]*exp(params[2]*(1/Tempdat[1]))
      icorr <- ((2*pi)/(3*lifedamparams[1]*pHdat[1]))*exp(-lifedamparams[2]/(K*Tempdat[1]))

      Dest <- ((r_0^3) + (Lifedat/(lifedamparams[1]*pHdat))*exp(-lifedamparams[3]/(K*Tempdat)))^(1/3)

      R2 <- 1 - sum(((Damdat^3) - (Dest^3))^2)/sum(((Damdat^3) - mean(Damdat^3))^2)
      return(list(lifedamparams,timepsuedo,icorr,R2))
    }
    damfit <- function(TimeDamfit,pH1,Temp1,params){
      ((r_0^3) + (TimeDamfit/(params[1]*pH1))*exp(-params[3]/(K*Temp1)))^(1/3)
    }
    table1labels <- c(Alabel1,"I_p0 (C/s)","E_a (eV)",time_output_names[1],icorrlabel1,"R\U00B2")
  }

  # Processing
  timefit <- linspace(0,max(data[,1]),1000)
  for(i in 1:length(unitnames)){
    Lifedam<-lifedamoutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],data[which(data[,3]==unitnames[i]),5],r_0)

    if(i==1){
      tableout<-c(Lifedam[[1]],Lifedam[[2]],Lifedam[[3]],Lifedam[[4]])
      damagefit<-damfit(timefit,stressvals[[2]][[i]],stressvals[[1]][[i]],Lifedam[[1]])
      timefit1<-timefit
    } else {
      tableout<-c(tableout,Lifedam[[1]],Lifedam[[2]],Lifedam[[3]],Lifedam[[4]])
      damagefit<-c(damagefit,NA,damfit(timefit,stressvals[[2]][[i]],stressvals[[1]][[i]],Lifedam[[1]]))
      timefit1<-c(timefit1,NA,timefit)
    }
  }
  tableout1<-matrix(tableout, nrow = length(unitnames), ncol = length(tableout)/length(unitnames), byrow = TRUE, dimnames = list(unitnames,table1labels))
  tableout2<-matrix(c(tableout1[,dim(tableout1)[2]-1],rep(1,length(unitnames)),unlist(stressvals)),nrow=length(stressvals[[1]]),ncol=2+stresscount,byrow=FALSE, dimnames = list(unitnames,time_output_names))

  # Return plot of degradation based on model
  # Data frame for data
  datastressname <- rep(0,dim(data)[1])
  for (i in 1:dim(data)[1]){
    datastressname[i] <- str_c(data[i,4]," ",names(data)[4],"/",data[i,5]," ",names(data)[5])
  }
  df1 <- data.frame(timedat = data[,1], pitdat = data[,2], group = datastressname)
  df2 <- data.frame(timefitvec = timefit1, damfitvec = damagefit)

  plotout<-ggplot() +
    geom_point(data=df1, aes(timedat,pitdat, colour=group, shape=group), size = 1.9) +
    geom_path(data=df2, aes(timefitvec,damfitvec), colour = "black", size = 0.5) +
    xlab(xlabel1) +
    ylab(ylabel1)

  return(list(corroutputtable = tableout1,stresstimetable = tableout2, corrplot = plotout))
}
