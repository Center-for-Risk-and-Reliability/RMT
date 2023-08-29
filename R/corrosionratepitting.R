# Pitting Corrosion Calculator
# Developed by Dr. Reuel Smith, 2022

corrosion.pitting <- function(data, Dam_0 = 0, Dam_cr = 1, matproperties = NULL, units=1) {
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

  # Case 1: Input data table (t, r, T) and material properties M, n, rho
  #         Output t vs. r plot, t vs. volume plot, I_p0 and E_a
  if (dim(data)[2] == 4 && length(matproperties) > 0){

  }

  # Case 2: Input data table (t, r, T, pH)
  #         Output t vs. r plot, t vs. volume plot, A and E_a
  if (dim(data)[2] == 5 && length(matproperties) == 0){
    # r^3 = r_0^3 + (Dt/(A x pH)) exp(-E_a/kT)
    # theta[1] ~ A, theta[2] ~ E_a
    lifedamoutput <- function(Lifedat,Damdat,Tempdat,pHdat,Dam_fail){
      params  <- pinv(matrix(c(rep(-1,length(Damdat)),-1/Tempdat),nrow = length(Damdat), ncol = 2, byrow = FALSE))%*%(log(Damdat^3) - log(Lifedat) + log(pHdat))
      lifedamparams <- c(exp(params[1]),params[2]*K)
      timepsuedo <- exp(params[1])*((Dam_cr^3) - (Dam_0^3))*pHdat[1]*exp(params[2]*(1/Tempdat[1]))

      Dest <- ((Dam_0^3) + (Lifedat/(lifedamparams[1]*pHdat))*exp(-lifedamparams[2]/(K*Tempdat)))^(1/3)

      R2 <- 1 - sum(((Damdat^3) - (Dest^3))^2)/sum(((Damdat^3) - mean(Damdat^3))^2)
      return(list(lifedamparams,timepsuedo,R2))
    }
    damfit <- function(TimeDamfit,pH1,Temp1,params){
      ((Dam_0^3) + (TimeDamfit/(params[1]*pH1))*exp(-params[2]/(K*Temp1)))^(1/3)
    }
  }

  # Case 3: Input data table (t, r, T, pH) and material properties M,n, and rho
  #         Output t vs. r plot, t vs. volume plot, I_p0, A and E_a
  if (dim(data)[2] == 5 && length(matproperties) > 0){

  }

  # Processing
  timefit <- linspace(0,max(data[,1]),1000)
  for(i in 1:length(unitnames)){
    Lifedam<-lifedamoutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],data[which(data[,3]==unitnames[i]),5],Dam_0)

    if(i==1){
      tableout<-c(Lifedam[[1]],Lifedam[[2]],Lifedam[[3]])
      damagefit<-damfit(timefit,stressvals[[2]][[i]],stressvals[[1]][[i]],Lifedam[[1]])
      timefit1<-timefit
    } else {
      tableout<-c(tableout,Lifedam[[1]],Lifedam[[2]],Lifedam[[3]])
      damagefit<-c(damagefit,NA,damfit(timefit,stressvals[[2]][[i]],stressvals[[1]][[i]],Lifedam[[1]]))
      timefit1<-c(timefit1,NA,timefit)
    }
  }
  tableout1<-matrix(tableout, nrow = length(unitnames), ncol = length(tableout)/length(unitnames), byrow = TRUE)
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
    # geom_point(data=df, aes(Nrunoff,Srunoff), colour = 'green4', shape=17, size = 1.9) +
    geom_path(data=df2, aes(timefitvec,damfitvec), colour = "black", size = 0.5) +
    # scale_x_continuous(trans = 'log10') +
    # scale_y_continuous(trans = 'log10') +
    # annotation_logticks() +
    xlab(names(data)[1]) +
    ylab(names(data)[2])

  return(list(tableout1,tableout2, plotout))
}
