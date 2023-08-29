# Least-Squares Accelerated Degradation Testing Estimator
# Developed by Dr. Reuel Smith, 2021-2022

adt.full.LSQ <- function(data,lifedam,D0,Tuse=293.15){
  # Load pracma library for pseudo-inverse
  library(pracma)
  library(dplyr)
  library(plotly)

  # Check the damage input.  If it is singular then treat it as such.
  # If it is a vector, make sure it is the same length as the column number
  # from the data (that is there is a unique damage condition for each unit)

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # Pulls names from column headers for the data
  col_names_data <- colnames(data)

  time_output_names <- c(col_names_data[1],"Censored",col_names_data[4:length(col_names_data)])

  if(lifedam=="Hamada" & missing(Tuse)){
    Tuse <- 293.15
  }


  # Check the damage input.  If it is singular then treat it as such.
  # If it is a vector, make sure it is the same length as the column number
  # from the data (that is there is a unique damage condition for each unit)

  if(length(D0)>1 && length(D0)<length(unitnames)){
    stop("'D0' has to be either a single value or a vector of the same length as the number of units in your data.")
  }

  # Also check if the damage is normalized or not.  If it is and the Mitsuom model is
  # entered, then it would run, but not if it is not.
  # if(max(data[,2])>=1 && lifedam=="Mitsuom"){
  #   stop("Select a different degradation model.  Mitsuom model is best reserved for damage between 0 and 1.")
  # }

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

  if(lifedam=="Linear"){
    # D = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(Damdat ~ poly(Lifedat, 1, raw=TRUE))
      lifedamparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      timepsuedo<-(Dam_fail - lifedamparams[1])/lifedamparams[2]
      R2 <- summary(params)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  if(lifedam=="Exponential"){
    # D = b*exp(a*t)
    # theta[1] ~ a, theta[2] ~ b
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(log(Damdat) ~ poly(Lifedat, 1, raw=TRUE))
      lifedamparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      timepsuedo<-(log(Dam_fail)-log(lifedamparams[2]))/lifedamparams[1]
      R2 <- summary(params)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  if(lifedam=="SquareRoot"){
    # D^(1/2) = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(sqrt(Damdat) ~ poly(Lifedat, 1, raw=TRUE))
      lifedamparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      timepsuedo<-(sqrt(Dam_fail)-lifedamparams[1])/lifedamparams[2]
      R2 <- summary(params)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  if(lifedam=="Power"){
    # D = b*(t^a)
    # theta[1] ~ a, theta[2] ~ b
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(log(Damdat) ~ poly(log(Lifedat), 1, raw=TRUE))
      lifedamparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      timepsuedo<-exp((log(Dam_fail) - log(lifedamparams[2]))/lifedamparams[1])
      R2 <- summary(params)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  if(lifedam=="Logarithmic"){
    # D = a + b*ln(t)
    # theta[1] ~ a, theta[2] ~ b
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(Damdat ~ poly(log(Lifedat), 1, raw=TRUE))
      lifedamparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      timepsuedo<-exp((Dam_fail - lifedamparams[1])/lifedamparams[2])
      R2 <- summary(params)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  # if(lifedam=="Gompertz"){
  #   # D = a + b^(c*t)
  #   # theta[1] ~ a, theta[2] ~ b, theta[3] ~ c
  #   lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
  #     a0<-min(Damdat) - 2
  #     b0c0<-pinv(matrix(rep(1,2*length(Damdat)),nrow = length(Damdat), ncol = 2))%*%(log(log(Damdat - a0))-log(Lifedat))
  #     b0<-exp(exp(b0c0[2]))
  #     c0<-exp(b0c0[1])
  #     datlabels<-colnames(data)
  #     params <- nls(datlabels[2] ~ a + b^(c*datlabels[1]), start = list(a = a0, b = b0, c = c0))
  #     lifedamparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
  #     timepsuedo<-exp((log(Dam_fail) - log(lifedamparams[2]))/lifedamparams[1])
  #     R2 <- summary(params)$r.squared
  #     return(list(lifedamparams,timepsuedo,R2))
  #   }
  # }

  if(lifedam=="LloydLipow"){
    # D = a - b/t
    # theta[1] ~ a, theta[2] ~ b
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(Damdat ~ poly((1/Lifedat), 1, raw=TRUE))
      lifedamparams <- c(summary(params)$coefficients[1,1],-summary(params)$coefficients[2,1])
      timepsuedo<-lifedamparams[2]/(lifedamparams[1] - Dam_fail)
      R2 <- summary(params)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  if(lifedam=="Mitsuom"){
    # D = 1/(1 + b*(t^a))
    # theta[1] ~ a, theta[2] ~ b
    # The Mitsuom model can't be solved by way of LSQ if there are degradation data greater than
    # or equal to 1, so the parameter analysis needs to be updated by way of MLE assuming
    # a lognormal fit (add a theta[3] for sigma_t).  A pre-fit of the data is set to an exponential
    # model though D(t) = EXP(m*t + b), then we use a lognormal updater
    # =========================================================================
    # Life damage output and update by MLE
    lifedamoutput <- function(Lifedat,Damdat,Dam_fail){
      # Compute LSQ estimate based on data not including Dam > 1
      params  <- lm(log(Damdat[which(Damdat<1)]^(-1) - 1) ~ poly(log(Lifedat[which(Damdat<1)]), 1, raw=TRUE))
      params_mb <- lm(log(Damdat) ~ poly(Lifedat, 1, raw=TRUE))
      lifedamparams0 <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      m_param<-summary(params_mb)$coefficients[2,1]
      b_param<-summary(params_mb)$coefficients[1,1]
      # lifedamparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      # Set the MSE equation
      loglik_Mitsuom <- function(theta){
        mean((Damdat - (1/(1 + theta[2]*(Lifedat^theta[1]))))^2)
      }
      lifedamparams_out <- nlm(loglik_Mitsuom, theta <- lifedamparams0, hessian=TRUE)
      lifedamparams <- lifedamparams_out$estimate[1:2]
      timepsuedo<-(((Dam_fail^(-1)) - 1)/lifedamparams[2])^(1/lifedamparams[1])
      Dammodel <- 1/(1 + lifedamparams[2]*(Lifedat^lifedamparams[1]))
      # R2 <- 1 - (sum((Damdat - Dammodel)^2)/sum((Damdat - mean(Dammodel))^2))
      R2 <- summary(params_mb)$r.squared
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  if(lifedam=="Hamada"){
    # D = 1/(1 + beta1*(t*exp(beta3*11605*(1/Tu - 1/Ti)))^beta2)
    # theta[1] ~ beta1, theta[2] ~ beta2, theta[3] ~ beta3
    lifedamoutput <- function(Lifedat,Damdat,Tempdat,Dam_fail){
      # Root out the data that exceeds D=1 so we fit the data properly
      iDam<-which(Damdat<1)
      params  <- pinv(matrix(c(rep(1,length(iDam)),log(Lifedat[iDam]),11605*((1/Tuse) - (1/Tempdat[iDam]))),nrow = length(iDam), ncol = 3, byrow = FALSE))%*%log((1/Damdat[iDam]) - 1)
      lifedamparams <- c(exp(params[1]),params[2],params[3]/params[2])
      timepsuedo<- exp((log((1/Dam_fail) - 1) - params[1] - params[3]*11605*((1/Tuse) - (1/Tempdat[iDam[1]])))/params[2])
      Dest<-1/(1 + lifedamparams[1]*((Lifedat[iDam]*exp(lifedamparams[3]*11605*((1/Tuse) - (1/Tempdat[iDam]))))^lifedamparams[2]))
      R2 <- 1 - sum((Damdat[iDam] - Dest)^2)/sum((Damdat[iDam] - mean(Damdat[iDam]))^2)
      return(list(lifedamparams,timepsuedo,R2))
    }
  }

  for(i in 1:length(unitnames)){
    if(length(D0)==1){
      if(lifedam=="Hamada"){
        Lifedam<-lifedamoutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],D0)
      } else {
        Lifedam<-lifedamoutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],D0)
      }
    } else if(length(D0)==length(unitnames)){
      if(lifedam=="Hamada"){
        Lifedam<-lifedamoutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],D0[i])
      } else {
        Lifedam<-lifedamoutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],D0[i])
      }
    }
    if(i==1){
      tableout<-c(Lifedam[[1]],Lifedam[[2]],Lifedam[[3]])
    } else {
      tableout<-c(tableout,Lifedam[[1]],Lifedam[[2]],Lifedam[[3]])
    }
  }
  tableout1<-matrix(tableout, nrow = length(unitnames), ncol = length(tableout)/length(unitnames), byrow = TRUE)
  if(stresscount==1){
    tableout2<-matrix(c(tableout1[,dim(tableout1)[2]-1],rep(1,length(stressvals)),stressvals),nrow=length(stressvals),ncol=3,byrow=FALSE, dimnames = list(unitnames,time_output_names))
  } else{
    tableout2<-matrix(c(tableout1[,dim(tableout1)[2]-1],rep(1,length(unitnames)),unlist(stressvals)),nrow=length(stressvals[[1]]),ncol=2+stresscount,byrow=FALSE, dimnames = list(unitnames,time_output_names))
  }
  return(list(tableout1,tableout2))
  # return(list(tableout1))
}
