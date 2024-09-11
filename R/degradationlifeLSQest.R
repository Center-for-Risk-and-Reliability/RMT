# Least-Squares Accelerated Degradation Testing Estimator
# Developed by Dr. Reuel Smith, 2021-2024

degradationlife.LSQest <- function(data,dl,dist="Lognormal",pp="Blom",D0,modelstress=NULL,xlabel=NULL,ylabel=NULL,Tuse=NULL){
  # Load pracma library for pseudo-inverse
  library(pracma)
  library(dplyr)
  library(plotly)
  library(stringr)
  library(ggplot2)

  # Legend colors
  col_legend <- c("red","blue","darkgreen","violet","aquamarine","orange","pink","darkblue","lightgreen","yellow","green")
  # Legend shapes
  shape_legend <- c(0:25)

  # UPDATE (2/12/2024) Now check to see if any time values (column 1) are zero with respect to
  # Power, Logarithmic, Lloyd-Lipow, and Mitsuom degradation life models (Add more where ln(L) is involved)
  data0 <- data
  if((dl=="Power" || dl=="Logarithmic" || dl=="LloydLipow" || dl=="Mitsuom") && min(data[,1])==0){
    # Nix time data at zero
    data <- data0[which(data0[,1]!=0),]
  }

  # UPDATE (2/6/2024 for you Mommy ♥♥♥)
  # Now we will take the data and reorganize it as if to obtain the LSQ estimates
  # for the degradation based on time and stress.  We will use this as a means to calculate
  # time-to-failure based on the degradation limit.
  # Sort the data first
  if(dim(data)[2] > 4){
    data_refit <- cbind(data[,2],rep(1,dim(data)[1]),data[,1],data[,4:dim(data)[2]])
  } else{
    data_refit <- cbind(data[,2],rep(1,dim(data)[1]),data[,1],data[,4])
  }

  # Set the column names as blank
  colnames(data_refit) <- rep(1,dim(data)[2])

  # Set names for the columns
  colnames(data_refit)[1] <- colnames(data)[2]
  colnames(data_refit)[2] <- "Censored"
  colnames(data_refit)[3] <- colnames(data)[1]
  if(dim(data)[2] > 4){
    colnames(data_refit)[4:dim(data)[2]] <- colnames(data)[4:dim(data)[2]]
  } else{
    colnames(data_refit)[4] <- colnames(data)[4]
  }
  xlabel1 <- colnames(data)[2]

  # Compute the LSQ estimates for the time/stress based degradation
  if (dist=="Weibull") {
    ppoutput <- probplot.wbl(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.wbl(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="Lognormal") {
    ppoutput <- probplot.logn(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.logn(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="Normal") {
    ppoutput <- probplot.nor(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.nor(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="Exponential") {
    ppoutput <- probplot.exp(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.exp(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="2PExponential") {
    ppoutput <- probplot.exp2P(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.exp2P(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="Gumbel") {
    ppoutput <- probplot.gumb(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.gumb(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="Logistic") {
    ppoutput <- probplot.logist(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.logist(data_refit,pp,xlabel1)$prob_plot
  }
  if (dist=="Loglogistic") {
    ppoutput <- probplot.loglogist(data_refit,pp,xlabel1)[[1]]
    plotoutput <- probplot.loglogist(data_refit,pp,xlabel1)$prob_plot
  }

  # return(ppoutput)

  # Check the damage input.  If it is singular then treat it as such.
  # If it is a vector, make sure it is the same length as the column number
  # from the data (that is there is a unique damage condition for each unit)

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # Pulls names from column headers for the data
  col_names_data <- colnames(data)

  # Pulls xlabel if given
  if(is.null(xlabel) == FALSE){
    col_names_data[1] <- xlabel
  }

  # Pulls ylabel if given
  if(is.null(ylabel) == FALSE){
    col_names_data[2] <- ylabel
  }

  time_output_names <- c(col_names_data[1],"Censored",col_names_data[4:length(col_names_data)])

  # NOTE: RCS02052024 Hold off on this for now until I reassign how the input will work
  # if(dl=="Hamada" & missing(Tuse)){
  #   Tuse <- 293.15
  # }

  # Check the damage input.  If it is singular then treat it as such.
  # If it is a vector, make sure it is the same length as the column number
  # from the data (that is there is a unique damage condition for each unit)

  if(length(D0)>1 && length(D0)<length(unitnames)){
    # NOTE: RCS02052024 Remember D0 can be a vector of length of the number of units as well as a single value
    stop("'D0' has to be either a single value or a vector of the same length as the number of units in your data.")
  }

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

  # Initialize parameter-stress relation if modelstress is not NULL (default)
  # First tick off which modelstress type it is based on the name (parameter-based - "1" or AF-based - "2")
  if(is.null(modelstress) == FALSE){
    if(modelstress == "Linear" || modelstress == "Exponential" || modelstress == "Exponential2" ||
       modelstress == "Arrhenius" || modelstress == "Eyring" || modelstress == "Eyring2" || modelstress == "Eyring3" || modelstress == "Eyring4" ||
       modelstress=="Power" || modelstress=="InversePower" || modelstress=="InversePower2" || modelstress=="Logarithmic" ||
       modelstress=="MultiStress" || modelstress=="TempHumidity" || modelstress=="TempNonthermal"){
      modelstresstype <- 1
    }
    if(modelstress == "LinearAF" || modelstress == "ExponentialAF" || modelstress == "ExponentialAF2" ||
       modelstress == "ArrheniusAF" || modelstress == "EyringAF" || modelstress == "EyringAF2" || modelstress == "EyringAF3" || modelstress == "EyringAF4"  ||
       modelstress=="PowerAF" || modelstress=="InversePowerAF" || modelstress=="InversePowerAF2" || modelstress=="LogarithmicAF" ||
       modelstress=="MultiStressAF" || modelstress=="TempHumidityAF" || modelstress=="TempNonthermalAF"){
      modelstresstype <- 2
    }
  } else{
    modelstresstype <- 0
  }

  if (is.null(modelstress) == FALSE && modelstress=="Linear"){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0

      params  <- lm(PARAM ~ poly(S, 1, raw=TRUE))
      psparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2] + params[1]*S
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"(b_0 + S*a_0)"
    logparam_txt<-"ln(b_0 + S*a_0)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Exponential"){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0, lsparams[3] - R^2
    paramstress <- function(PARAM,S){
      # theta[1] ~ a, theta[2] ~ b
      # If stressmodel is based on parameter, it is based on b

      params  <- lm(log(PARAM) ~ poly(S, 1, raw=TRUE))
      psparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2]*exp(params[1]*S)
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"b_0*exp(a_0*S)"
    logparam_txt<-"(log(b_0) + a_0*S)"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Exponential2"){
    # lsparams[1] - parameter a_0, lsparams[2] - parameter b_0
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0

      params  <- lm(log(PARAM) ~ poly(1/S, 1, raw=TRUE))
      psparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2]*exp(params[1]/S)
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"b_0*exp(a_0/S)"
    logparam_txt<-"(log(b_0) + a_0/S)"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Arrhenius"){
    # psparams[1] - parameter Ea, psparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    paramstress <- function(PARAM,S){
      # theta[1] ~ E_a_0, theta[2] ~ b_0
      params  <- lm(log(PARAM) ~ poly(1/S, 1, raw=TRUE))
      psparams <- c(K*summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2]*exp(params[1]/(K*S))
    }
    # Writeup for the output text
    params_txt<-c("E_a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"b_0*exp(E_a_0)/(K*S))"
    logparam_txt<-"(log(b_0) + (E_a_0/(K*S)))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Eyring"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0
      Lvals<-log(PARAM) + log(S)
      Svals<-matrix(c(1/S,rep(1,length(S))),nrow=length(S),ncol=2,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      psparams <- c(params)
      psparams[2]<-exp(psparams[2])

      return(psparams)
    }
    stressparamfit <- function(params,S){
      (params[2]/S)*exp(params[1]/S)
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"(b_0/S)*exp(a_0/S)"
    logparam_txt<-"(log(b_0) - log(S) + (a_0/S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Eyring2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0
      Lvals<-log(PARAM) + log(S)
      Svals<-matrix(c(rep(-1,length(S)), 1/S),nrow=length(S),ncol=2,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      psparams <- c(params)

      return(psparams)
    }
    stressparamfit <- function(params,S){
      (1/S)*exp(-(params[1] - (params[2]/S)))
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-"Eyring (Type-2)"
    param_txt2<-"(1/S)*exp(-(a_0 - (b_0/S)))"
    logparam_txt<-"(-log(S) - a_0 + (b_0/S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Power"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0
      params  <- lm(log(PARAM) ~ poly(log(S), 1, raw=TRUE))
      psparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2]*(S^params[1])
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"b_0*(S^a_0)"
    logparam_txt<-"(ln(b_0) + a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="InversePower"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0
      params  <- lm(log(PARAM) ~ poly(log(S), 1, raw=TRUE))
      psparams <- c(-summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2]*(S^-params[1])
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-"Inverse Power"
    param_txt2<-"b_0*(S^-a_0)"
    logparam_txt<-"(ln(b_0) - a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="InversePower2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0
      params  <- lm(log(PARAM) ~ poly(log(S), 1, raw=TRUE))
      psparams <- c(-summary(params)$coefficients[2,1],exp(-summary(params)$coefficients[1,1]))

      return(psparams)
    }
    stressparamfit <- function(params,S){
      1/(params[2]*(S^params[1]))
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-"Inverse Power"
    param_txt2<-"1/[b_0*(S^a_0)]"
    logparam_txt<-"(-ln(b_0) - a_0ln(S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Logarithmic"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0
      params  <- lm(PARAM ~ poly(log(S), 1, raw=TRUE))
      psparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[2] + params[1]*log(S)
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0")
    ps_txt<-modelstress
    param_txt2<-"(b_0 + a_0*ln(S))"
    logparam_txt<-"ln(b_0 + a_0*ln(S))"
  }
  if (is.null(modelstress) == FALSE && modelstress=="MultiStress"){
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),S),nrow=length(ppoutput)/3,ncol=1+length(ppoutput[[1]]),byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    lnLmodel <- Svals%*%lsparams
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    # Writeup for the output text
    params_txt<-paste("a_",c(0:length(S[1,])),sep="")
    ps_txt<-"Multi-Stress"
    param_txt2<-"exp(a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n)"
    logparam_txt<-"a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n"
  }
  if (is.null(modelstress) == FALSE && modelstress=="TempHumidity"){
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    if((dim(data)[2]-3)<2) {
      stop('Select a data set with more than one stress type.')
    }
    paramstress <- function(PARAM,S){
      # theta[1] ~ E_a_0, theta[2] ~ b_0
      Lvals<-log(PARAM)
      Svals<-matrix(c(rep(1,length(S[,1])),1/S[,1],1/S[,2]),nrow=length(unitnames),ncol=3,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      params[1]<-exp(params[1])
      psparams <- c(params)

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[1]*exp((params[2]/S[,1]) + (params[3]/S[,2]))
    }
    # Writeup for the output text
    params_txt<-c("A_0","a_0","b_0")
    ps_txt<-"Temperature-Humidity"
    param_txt2<-"A_0 exp(a_0/S + b_0/H)"
    logparam_txt<-"ln(A_0) + a_0/S + b_0/H"
  }

  if (is.null(modelstress) == FALSE && modelstress=="TempNonthermal"){
    if((dim(data)[2]-3)<2) {
      stop('Select a data set with more than one stress type.')
    }
    paramstress <- function(PARAM,S){
      # theta[1] ~ a_0, theta[2] ~ b_0, theta[3] ~ c_0
      Lvals<-log(PARAM)
      Svals<-matrix(c(1/S[,1],-log(S[,2]),rep(1,length(S[,1]))),nrow=length(unitnames),ncol=3,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      params[3]<-exp(params[3])
      psparams <- c(params)

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[3]/((S[2]^params[2])*exp(-params[1]/S[1]))
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0","c_0")
    ps_txt<-"Temperature-Non-thermal"
    param_txt2<-"c_0/(U^b_0 * exp(-a_0/S))"
    logparam_txt<-"a_0(1/S) - b_0*ln(U) + ln(c_0)"
  }

  if (is.null(modelstress) == FALSE && modelstress=="Eyring3"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    if((dim(data)[2]-3)<2) {
      stop('Select a data set with more than one stress type.')
    }
    paramstress <- function(PARAM,S){
      # lsparams[1] - parameter a, lsparams[2] - parameter b
      # lsparams[3] - parameter c, lsparams[4] - parameter d
      Lvals<-log(PARAM)+log(S[,1])
      Svals<-matrix(c(rep(1,length(S[,1])),1/S[,1],S[,2],S[,2]/S[,1]),nrow=length(unitnames),ncol=4,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      psparams <- c(params)

      return(psparams)
    }
    stressparamfit <- function(params,S){
      (1/S[1])*exp((params[1] + (params[2]/S[1])) + ((params[3] + (params[4]/S[1]))*S[2]))
    }
    # Writeup for the output text
    params_txt<-c("a_0","b_0","c_0","d_0")
    ps_txt<-"Eyring (Type 3)"
    param_txt2<-"(1/S) exp((a_0 + (b_0/S)) + (c_0 + (d_0/S)) U)"
    logparam_txt<-"-ln(S) + (a_0 + (b_0/S)) + (c_0 + (d_0/S)) U"
  }
  if (is.null(modelstress) == FALSE && modelstress=="Eyring4"){
    # lsparams[1] - parameter A, lsparams[2] - parameter b
    # lsparams[3] - parameter Ea
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    if((dim(data)[2]-3)<2) {
      stop('Select a data set with more than one stress type.')
    }
    paramstress <- function(PARAM,S){
      # lsparams[1] - parameter A, lsparams[2] - parameter b
      # lsparams[3] - parameter Ea
      Lvals<-log(PARAM)
      Svals<-matrix(c(rep(1,length(S[,1])),-log(S[,2]),1/S[,1]),nrow=length(unitnames),ncol=3,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      psparams <- c(params)
      psparams[1]<-exp(psparams[1])
      psparams[3]<-K*psparams[3]

      return(psparams)
    }
    stressparamfit <- function(params,S){
      params[1]*exp(params[3]/(K*S[1]))*(S[2]^-params[2])
    }

    params_txt<-c("A_0","b_0","E_a_0")
    ps_txt<-"Eyring (Type 3)"
    param_txt2<-"A_0 exp(E_a_0/(K*S)) U^-b_0"
    logparam_txt<-"ln(A_0) + (E_a_0/(K*S)) - b_0 ln(U)"
  }

  # Initialize AF-stress relation if modelstress is not NULL
  # INitial test will be for the TempNontermalAF model for several scenarios but mostly for Power

  # Degradation-life model setups.  Least-squares fit for standard degradation life models
  # without stress correlation.
  if(dl=="Linear"){
    # D = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(Damdat ~ poly(Lifedat, 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      timepsuedo<-(Dam_fail - dlparams[1])/dlparams[2]
      R2 <- summary(params)$r.squared
      SSE <- sum((fitted(params) - Damdat)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      params[1] + TimeDamfit*params[2]
    }
    # Writeup for the output text
    dl_txt <- dl
    dlmodel_txt0<-"a + bl"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c(params_txt,"b")
        paramref_txt<-"a(S) = "
        paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
        params_txt2<-c("\U03BC_b","\U03C3_b")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  if(dl=="Exponential"){
    # D = b*exp(a*t)
    # theta[1] ~ a, theta[2] ~ b
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(log(Damdat) ~ poly(Lifedat, 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      timepsuedo<-(log(Dam_fail)-log(dlparams[2]))/dlparams[1]
      R2 <- summary(params)$r.squared
      SSE <- sum((exp(fitted(params)) - Damdat)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      params[2]*exp(TimeDamfit*params[1])
    }
    # Writeup for the output text
    dl_txt <- dl
    dlmodel_txt0<-"b exp(al)"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c("a",params_txt)
        paramref_txt<-"b(S) = "
        paramref_txt2<-"a ~ NOR(\U03BC_a,\U03C3_a)"
        params_txt2<-c("\U03BC_a","\U03C3_a")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  if(dl=="SquareRoot"){
    # D^(1/2) = a + b*t
    # theta[1] ~ a, theta[2] ~ b
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(sqrt(Damdat) ~ poly(Lifedat, 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      timepsuedo<-(sqrt(Dam_fail)-dlparams[1])/dlparams[2]
      R2 <- summary(params)$r.squared
      SSE <- sum(((fitted(params)^2) - Damdat)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      (params[1] + TimeDamfit*params[2])^2
    }
    # Writeup for the output text
    dl_txt <- "Square-Root"
    dlmodel_txt0<-"(a + bl)\U00B2"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c(params_txt,"b")
        paramref_txt<-"a(S) = "
        paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
        params_txt2<-c("\U03BC_b","\U03C3_b")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  # Need to work in a condition here for when degradation limit is more for %
  if(dl=="Power"){
    # D = b*(t^a)
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      # theta[1] ~ a, theta[2] ~ b
      # If stressmodel is based on parameter, it is based on b

      params  <- lm(log(Damdat) ~ poly(log(Lifedat), 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))

      timepsuedo<-exp((log(Dam_fail) - log(dlparams[2]))/dlparams[1])
      R2 <- summary(params)$r.squared
      SSE <- sum((exp(fitted(params)) - Damdat)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      params[2]*(TimeDamfit^params[1])
    }

    damfit1 <- function(TimeDamfit,params,S){
      stressparamfit(params[2:length(params)],S)*(TimeDamfit^params[1])
    }
    # if (is.null(modelstress) == FALSE && modelstress=="TempNonthermalAF"){
    #   dloutput2 <- function(Lifedat,Damdat,Dam_fail){
    #     # theta[1] ~ a, theta[2] ~ b, theta[3] ~ a_0, theta[4] ~ b_0
    #
    #     Lvals<-log(Damdat)
    #     Svals<-matrix(c(log(Lifedat),rep(1,length(S[,1])),log(S[,2]/Tuse[2]),((1/Tuse[1]) - (1/S[,1]))),nrow=length(unitnames),ncol=3,byrow=FALSE)
    #     params  <- pinv(Svals)%*%Lvals
    #     dlparams <- c(params)
    #     dlparams[2]<-exp(dlparams[2])
    #     dlparams[3]<-dlparams[4]/dlparams[1]
    #     dlparams[4]<-dlparams[3]/dlparams[1]
    #
    #     params  <- lm(log(Damdat) ~ poly(log(Lifedat), 1, raw=TRUE))
    #     dlparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    #
    #     timepsuedo<-exp((log(Dam_fail) - log(dlparams[2]))/dlparams[1])
    #     R2 <- summary(params)$r.squared
    #     return(list(dlparams,timepsuedo,R2))
    #   }
    # }
    # Writeup for the output text
    dl_txt <- dl
    dlmodel_txt0<-"bl\U00AA"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c("a",params_txt)
        paramref_txt<-"b(S) = "
        paramref_txt2<-"a ~ NOR(\U03BC_a,\U03C3_a)"
        params_txt2<-c("\U03BC_a","\U03C3_a")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  if(dl=="Logarithmic"){
    # D = a + b*ln(t)
    # theta[1] ~ a, theta[2] ~ b
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(Damdat ~ poly(log(Lifedat), 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      timepsuedo<-exp((Dam_fail - dlparams[1])/dlparams[2])
      R2 <- summary(params)$r.squared
      SSE <- sum((fitted(params) - Damdat)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      params[1] + params[2]*log(TimeDamfit)
    }
    damfit1 <- function(TimeDamfit,params,S){
      stressparamfit(params[1:(length(params)-1)],S) + params[length(params)]*log(TimeDamfit)
    }
    # Writeup for the output text
    dl_txt <- dl
    dlmodel_txt0<-"a + b ln(l)"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c(params_txt,"b")
        paramref_txt<-"a(S) = "
        paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
        params_txt2<-c("\U03BC_b","\U03C3_b")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  # if(dl=="Gompertz"){
  #   # D = a + b^(c*t)
  #   # theta[1] ~ a, theta[2] ~ b, theta[3] ~ c
  #   dloutput <- function(Lifedat,Damdat,Dam_fail){
  #     a0<-min(Damdat) - 2
  #     b0c0<-pinv(matrix(rep(1,2*length(Damdat)),nrow = length(Damdat), ncol = 2))%*%(log(log(Damdat - a0))-log(Lifedat))
  #     b0<-exp(exp(b0c0[2]))
  #     c0<-exp(b0c0[1])
  #     datlabels<-colnames(data)
  #     params <- nls(datlabels[2] ~ a + b^(c*datlabels[1]), start = list(a = a0, b = b0, c = c0))
  #     dlparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
  #     timepsuedo<-exp((log(Dam_fail) - log(dlparams[2]))/dlparams[1])
  #     R2 <- summary(params)$r.squared
  #     return(list(dlparams,timepsuedo,R2))
  #   }
  # }

  if(dl=="LloydLipow"){
    # D = a - b/t
    # theta[1] ~ a, theta[2] ~ b
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(Damdat ~ poly((1/Lifedat), 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[1,1],-summary(params)$coefficients[2,1])
      timepsuedo<-dlparams[2]/(dlparams[1] - Dam_fail)
      R2 <- summary(params)$r.squared
      SSE <- sum((fitted(params) - Damdat)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      params[1] - (params[2]/TimeDamfit)
    }
    # Writeup for the output text
    dl_txt <- "Lloyd-Lipow"
    dlmodel_txt0<-"a - b/l"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c(params_txt,"b")
        paramref_txt<-"a(S) = "
        paramref_txt2<-"b ~ NOR(\U03BC_b,\U03C3_b)"
        params_txt2<-c("\U03BC_b","\U03C3_b")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  if(dl=="Mitsuom"){
    # D = 1/(1 + b*(t^a))
    # theta[1] ~ a, theta[2] ~ b
    # The Mitsuom model can't be solved by way of LSQ if there are degradation data greater than
    # or equal to 1, so the parameter analysis needs to be updated by way of MLE assuming
    # a lognormal fit (add a theta[3] for sigma_t).  A pre-fit of the data is set to an exponential
    # model though D(t) = EXP(m*t + b), then we use a lognormal updater
    # =========================================================================
    # Life damage output and update by MLE
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      # Compute LSQ estimate based on data not including Dam > 1
      params  <- lm(log(Damdat[which(Damdat<1)]^(-1) - 1) ~ poly(log(Lifedat[which(Damdat<1)]), 1, raw=TRUE))
      params_mb <- lm(log(Damdat) ~ poly(Lifedat, 1, raw=TRUE))
      dlparams0 <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      m_param<-summary(params_mb)$coefficients[2,1]
      b_param<-summary(params_mb)$coefficients[1,1]
      # dlparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      # Set the MSE equation
      loglik_Mitsuom <- function(theta){
        mean((Damdat - (1/(1 + theta[2]*(Lifedat^theta[1]))))^2)
      }
      dlparams_out <- nlm(loglik_Mitsuom, theta <- dlparams0, hessian=TRUE)
      dlparams <- dlparams_out$estimate[1:2]
      timepsuedo<-(((Dam_fail^(-1)) - 1)/dlparams[2])^(1/dlparams[1])
      Dammodel <- 1/(1 + dlparams[2]*(Lifedat^dlparams[1]))
      # R2 <- 1 - (sum((Damdat - Dammodel)^2)/sum((Damdat - mean(Dammodel))^2))
      R2 <- summary(params_mb)$r.squared
      SSE <- sum((Damdat - Dammodel)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,params){
      1/(1 + params[2]*(TimeDamfit^params[1]))
    }
    # Writeup for the output text
    dl_txt <- dl
    dlmodel_txt0<-"(1 + bl\U00AA)\U207B\U00B9"
    if(is.null(modelstress) == FALSE){
      if(modelstresstype == 1){
        dlparams_txt<-c("a","b")
        # dlparams_txt<-c("a",params_txt)
        paramref_txt<-"b(S) = "
        paramref_txt2<-"a ~ NOR(\U03BC_a,\U03C3_a)"
        params_txt2<-c("\U03BC_a","\U03C3_a")
      }
      if(modelstresstype == 2){
        dlparams_txt<-c("a","b")
      }
    } else{
      dlparams_txt<-c("a","b")
    }
  }

  if(dl=="Hamada"){
    # D = 1/(1 + beta1*(t*exp(beta3*11605*(1/Tu - 1/Ti)))^beta2)
    # theta[1] ~ beta1, theta[2] ~ beta2, theta[3] ~ beta3
    dloutput <- function(Lifedat,Damdat,Tempdat,Dam_fail){
      # Root out the data that exceeds D=1 so we fit the data properly
      iDam<-which(Damdat<1)
      params  <- pinv(matrix(c(rep(1,length(iDam)),log(Lifedat[iDam]),11605*((1/Tuse) - (1/Tempdat[iDam]))),nrow = length(iDam), ncol = 3, byrow = FALSE))%*%log((1/Damdat[iDam]) - 1)
      dlparams <- c(exp(params[1]),params[2],params[3]/params[2])
      timepsuedo<- exp((log((1/Dam_fail) - 1) - params[1] - params[3]*11605*((1/Tuse) - (1/Tempdat[iDam[1]])))/params[2])
      Dest<-1/(1 + dlparams[1]*((Lifedat[iDam]*exp(dlparams[3]*11605*((1/Tuse) - (1/Tempdat[iDam]))))^dlparams[2]))
      R2 <- 1 - sum((Damdat[iDam] - Dest)^2)/sum((Damdat[iDam] - mean(Damdat[iDam]))^2)
      SSE <- sum((Damdat[iDam] - Dest)^2)
      return(list(dlparams,timepsuedo,R2,SSE))
    }
    damfit <- function(TimeDamfit,TempDamfit,params){
      1/(1 + params[1]*(TimeDamfit*exp(params[3]*11605*(1/Tuse - 1/TempDamfit)))^params[2])
    }
    # Writeup for the output text
    dl_txt <- dl
    dlmodel_txt0<-"(1 + bl\U00AA)\U207B\U00B9"
    dlparams_txt<-c("a","b")
  }

  # life damage model "CrackProp1" is an assumed effective zero initial crack length
  # model adjusted for an initial crack length of 1 mm or 0.001 m.  It was used in the
  # ENRE 641 2022 Final Exam with an assumed DS of 200 MPa.
  # NOTE: Prepare data first by converting to meters and then subtracting 0.001 to zero it
  if(dl=="CrackProp1"){
    # D = a_0(=0.001 m) + exp((2/(2-m))*LN(N) + (2/(2-m))*(LN(C*sqrt-pi) + LN(1 - 0.5*m) + LN(DS)))
    # theta[1] ~ C, theta[2] ~ m
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(log(Damdat) ~ poly(log(Lifedat), 1, raw=TRUE))
      gamparam <- 1/summary(params)$coefficients[2,1]
      Aparam <- exp((summary(params)$coefficients[1,1]/summary(params)$coefficients[2,1]) - log(gamparam) - log(200))
      dlparams <- c(Aparam/sqrt(pi),2*(1-gamparam))
      timepsuedo<-((Dam_fail - 0.001)^(1 - 0.5*dlparams[2]))/(dlparams[1]*(1 - 0.5*dlparams[2])*200*sqrt(pi))
      R2 <- summary(params)$r.squared
      return(list(dlparams,timepsuedo,R2))
    }
    damfit <- function(TimeDamfit,params){
      0.001 + exp((2/(2 - params[2]))*log(TimeDamfit) + (2/(2 - params[2]))*(log(params[1]*sqrt(pi)) + log(1 - 0.5*params[2]) + log(200)))
    }
  }
  if(dl=="CrackProp2"){
    # D = a_0(=0.001 m)*exp(pi*C*DS^2*N_a)
    # theta[1] ~ logC (C is too small to control so take the natural log and process that)
    dloutput <- function(Lifedat,Damdat,Dam_fail){
      params  <- lm(log(Damdat/0.001) ~ 0 + poly(Lifedat, 1, raw=TRUE))
      dlparams <- c(summary(params)$coefficients[1,1]/(pi*(200^2)))
      timepsuedo<-(1/(pi*dlparams*(200^2)))*log(Dam_fail/0.001)
      R2 <- summary(params)$r.squared
      return(list(dlparams,timepsuedo,R2))
    }
    damfit <- function(TimeDamfit,params){
      0.001*exp(pi*params[1]*(200^2)*TimeDamfit)
    }
  }
  # if(dl=="KondoWei"){
  #   # r^3 = r_0^3 + (Dt/(A x pH)) exp(-E_a/kT)
  #   # theta[1] ~ A, theta[2] ~ E_a
  #
  #   # dloutput <- function(Lifedat,Damdat,Dam_fail){
  #   #   params  <- lm(log(Damdat) ~ poly(log(Lifedat), 1, raw=TRUE))
  #   #   dlparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
  #   #   timepsuedo<-exp((log(Dam_fail) - log(dlparams[2]))/dlparams[1])
  #   #   R2 <- summary(params)$r.squared
  #   #   return(list(dlparams,timepsuedo,R2))
  #   # }
  #   # damfit <- function(TimeDamfit,params){
  #   #   params[2]*(TimeDamfit^params[1])
  #   # }
  #   dloutput <- function(Lifedat,Damdat,Tempdat,pHdat,Dam_fail){
  #     params  <- pinv(matrix(c(rep(-1,length(Damdat)),-1/Tempdat),nrow = length(Damdat), ncol = 2, byrow = FALSE))%*%(log(Damdat^3) - log(Lifedat) + log(pHdat))
  #     dlparams <- c(exp(params[1]),params[2]*K)
  #     timepsuedo <- exp(params[1])*((Dam_cr^3) - (Dam_0^3))*pHdat[1]*exp(params[2]*(1/Tempdat[1]))
  #     Dest <- ((Dam_0^3) + (Lifedat/(dlparams[1]*pHdat))*exp(-dlparams[2]/(K*Tempdat)))^(1/3)
  #     R2 <- 1 - sum(((Damdat^3) - (Dest^3))^2)/sum(((Damdat^3) - mean(Damdat^3))^2)
  #     return(list(dlparams,timepsuedo,R2))
  #   }
  #   damfit <- function(TimeDamfit,pH1,Temp1,params){
  #     ((Dam_0^3) + (TimeDamfit/(params[1]*pH1))*exp(-params[2]/(K*Temp1)))^(1/3)
  #   }
  # }


  # Processing

  # Generate model parameters and  best-fit curves by unit

  for(i in 1:length(unitnames)){
    if(length(D0)==1){
      if(dl=="Hamada"){
        Lifedam<-dloutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],D0)
      } else if(dl=="KondoWei"){
        Lifedam<-dloutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],data[which(data[,3]==unitnames[i]),5],D0)
      } else{
        Lifedam<-dloutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],D0)
      }
    } else if(length(D0)==length(unitnames)){
      if(dl=="Hamada"){
        Lifedam<-dloutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],data[which(data[,3]==unitnames[i]),4],D0[i])
      } else {
        Lifedam<-dloutput(data[which(data[,3]==unitnames[i]),1],data[which(data[,3]==unitnames[i]),2],D0[i])
      }
    }
    if(i==1){
      tableout<-c(Lifedam[[1]],Lifedam[[2]],Lifedam[[3]],Lifedam[[4]])
    } else{
      tableout<-c(tableout,Lifedam[[1]],Lifedam[[2]],Lifedam[[3]],Lifedam[[4]])
    }
  }
  tableout1<-matrix(tableout, nrow = length(unitnames), ncol = length(tableout)/length(unitnames), byrow = TRUE)

  # return(list(modelstresstype,ppoutput,tableout1,stressvals))
  tableout0 <- tableout1

  # If parameter based stress model then recompute based on stress and parameter (RCS02082024)
  if(modelstresstype == 1 && (dl=="Linear" || dl=="SquareRoot" || dl=="Logarithmic" || dl=="LloydLipow")){
    # Parameter fit to "a"
    if((dim(data)[2]-3)<2) {
      params_0 <- paramstress(tableout1[,1],stressvals)
    } else{
      params_0 <- paramstress(tableout1[,1],matrix(unlist(stressvals),nrow = length(unitnames), ncol = (dim(data)[2]-3), byrow = TRUE))
    }
    params_1 <- c(probplot.nor(cbind(tableout1[,2],rep(1,length(unitnames)),rep(1,length(unitnames))),pp="Blom")[[1]][[2]])
  }
  if(modelstresstype == 1 && (dl=="Exponential" || dl=="Power" || dl=="Mitsuom")){
    # Parameter fit to "b"
    if((dim(data)[2]-3)<2) {
      params_0 <- paramstress(tableout1[,2],stressvals)
    } else{
      params_0 <- paramstress(tableout1[,2],matrix(unlist(stressvals),nrow = length(unitnames), ncol = (dim(data)[2]-3), byrow = TRUE))
    }
    params_1 <- c(probplot.nor(cbind(tableout1[,1],rep(1,length(unitnames)),rep(1,length(unitnames))),pp="Blom")[[1]][[2]])
  }

  # return(list(modelstresstype,tableout1,dlparams_txt,params_0))

  # Return plot of degradation based on model
  # Data frame for data
  datastressname <- rep(0,dim(data0)[1])
  if((dim(data)[2]-3)<2){
    for (i in 1:dim(data0)[1]){
      datastressname[i] <- str_c(data0[i,4]," ",names(data0)[4])
    }
  } else{
    for (i in 1:dim(data0)[1]){
      datastressname[i] <- str_c(data0[i,4]," ",names(data0)[4],"/",data0[i,5]," ",names(data0)[5])
    }
  }

  # return(tableout1[,3][1])

  # Generate the best-fit curves by unit (Now separate operations)
  for(i in 1:length(unitnames)){
    # Establish time fit vector for the final plot up to the maximum time value in the data
    # 5/6/2024 - Time fit now adheres to the endurance limit as its stopping point
    if(dl=="Power" || dl=="Logarithmic" || dl=="LloydLipow" || dl=="CrackProp1"){
      timefit <- linspace(0.01,tableout1[,3][i],1000)
    } else{
      timefit <- linspace(0,tableout1[,3][i],1000)
    }

    if(i==1){
      if(dl=="Hamada"){
        damagefit<-damfit(timefit,rep(data[which(data[,3]==unitnames[i]),4][1],1000),tableout1[i,1:2])
        timefit1<-timefit
        datastressname2<-rep(datastressname[which(data0[,3]==unitnames[i])[1]],1000)
      } else if(dl=="KondoWei"){
        damagefit<-damfit(timefit,rep(data[which(data[,3]==unitnames[i]),4][1],1000),rep(data[which(data[,3]==unitnames[i]),5][1],1000),tableout1[i,1:2])
        timefit1<-timefit
        datastressname2<-rep(datastressname[which(data0[,3]==unitnames[i])[1]],1000)
      } else {
        if(modelstresstype == 0){
          damagefit<-damfit(timefit,tableout1[i,1:2])
        }
        if(modelstresstype == 1 && is.list(stressvals) == FALSE){
          damagefit<-damfit(timefit,tableout1[i,1:2])
          # damagefit<-damfit1(timefit,tableout1[i,1:length(dlparams_txt)],stressvals[i])
        }
        if(modelstresstype == 1 && is.list(stressvals) == TRUE){
          damagefit<-damfit(timefit,tableout1[i,1:2])
          # damagefit<-damfit1(timefit,tableout1[i,1:length(dlparams_txt)],c(stressvals[[1]][[i]],stressvals[[2]][[i]]))
        }
        timefit1<-timefit
        datastressname2<-rep(datastressname[which(data0[,3]==unitnames[i])[1]],1000)
      }
    } else {
      if(dl=="Hamada"){
        if(modelstresstype == 0){
          damagefit<-c(damagefit,NA,damfit(timefit,rep(data[which(data[,3]==unitnames[i]),4][1],1000),tableout1[i,1:2]))
        }
        timefit1<-c(timefit1,NA,timefit)
        datastressname2<-c(datastressname2,datastressname2[length(datastressname2)],rep(datastressname[which(data0[,3]==unitnames[i])[1]],1000))
      } else{
        if(modelstresstype == 0){
          damagefit<-c(damagefit,NA,damfit(timefit,tableout1[i,1:2]))
        }
        if(modelstresstype == 1 && is.list(stressvals) == FALSE){
          damagefit<-c(damagefit,NA,damfit(timefit,tableout1[i,1:2]))
          # damagefit<-c(damagefit,NA,damfit1(timefit,tableout1[i,1:length(dlparams_txt)],stressvals[i]))
        }
        if(modelstresstype == 1 && is.list(stressvals) == TRUE){
          damagefit<-c(damagefit,NA,damfit(timefit,tableout1[i,1:2]))
          # damagefit<-c(damagefit,NA,damfit1(timefit,tableout1[i,1:length(dlparams_txt)],c(stressvals[[1]][[i]],stressvals[[2]][[i]])))
        }
        timefit1<-c(timefit1,NA,timefit)
        datastressname2<-c(datastressname2,datastressname2[length(datastressname2)],rep(datastressname[which(data0[,3]==unitnames[i])[1]],1000))
      }
    }
  }

  # return(tableout1)

  # Forms first table containing the output parameters and coefficient of determination by unit
  if(stresscount==1){
    tableout2<-matrix(c(tableout1[,dim(tableout1)[2]-2],rep(1,length(stressvals)),stressvals),nrow=length(stressvals),ncol=3,byrow=FALSE, dimnames = list(unitnames,time_output_names))
  } else{
    tableout2<-matrix(c(tableout1[,dim(tableout1)[2]-2],rep(1,length(unitnames)),unlist(stressvals)),nrow=length(stressvals[[1]]),ncol=2+stresscount,byrow=FALSE, dimnames = list(unitnames,time_output_names))
  }

  if(dl=="CrackProp1"){
    df1 <- data.frame(timedat = data0[,1], damdat = data0[,2]+0.001, group = datastressname)
  } else{
    df1 <- data.frame(timedat = data0[,1], damdat = data0[,2], group = datastressname)
  }

  df2 <- data.frame(timefitvec = timefit1, damfitvec = damagefit, group = datastressname2)
  df3 <- data.frame(timefitvec = timefit1, damfitvec = rep(D0,length(timefit1)), endurance = rep("Degradation Limit",length(timefit1)))
  # df2 <- data.frame(timefitvec = timefit1, damfitvec = damagefit, group = datastressname2)

  # return(df2)

  plotout<-ggplot() +
    geom_point(data=df1, aes(timedat,damdat, colour=group, shape=group), size = 2.9) +
    # geom_point(data=df, aes(Nrunoff,Srunoff), colour = 'green4', shape=17, size = 1.9) +
    geom_path(data=df2, aes(timefitvec,damfitvec, linetype = group), colour = "black", size = 0.5) +
    scale_shape_manual(values=shape_legend[1:length(unitnames)]) +
    scale_color_manual(values=col_legend[1:length(unitnames)])+
    # scale_x_continuous(trans = 'log10') +
    # scale_y_continuous(trans = 'log10') +
    # annotation_logticks() +
    xlab(col_names_data[1]) +
    ylab(col_names_data[2])
  plotout<-plotout + geom_path(data=df3, aes(timefitvec,damfitvec), linetype = "dashed", colour = "red", size = 0.9)

  # Produce some output text that summarizes the results
  if(modelstresstype == 0){
    cat(c("Least-Squares estimates for the ",dl_txt," Degradation-Life model.\n\nD(l) = ",dlmodel_txt0,"\n\n"),sep = "")
    print(matrix(tableout, nrow = length(unitnames), ncol = length(tableout)/length(unitnames), byrow = TRUE,dimnames = list(unitnames,c(dlparams_txt,"Pseudo-Failure Times","R\U00B2","SSE"))))
    cat("\n")

    return(list(tableout1,tableout2,plotout))
  }
  # return(list(dl_txt,ps_txt,dlmodel_txt0,paramref_txt,param_txt2,paramref_txt2))
  if(modelstresstype == 1 || modelstresstype == 2){
    cat(c("Least-Squares estimates for the ",dl_txt,"-",ps_txt," Degradation-Life-Stress model.\n\nD(l,S) = ",dlmodel_txt0," where ",paramref_txt,param_txt2," and ",paramref_txt2,"\n\n"),sep = "")
    print(matrix(tableout, nrow = length(unitnames), ncol = length(tableout)/length(unitnames), byrow = TRUE,dimnames = list(unitnames,c(dlparams_txt,"Pseudo-Failure Times","R\U00B2","SSE"))))
    cat("\n")
    print(matrix(c(params_0,params_1), nrow = length(params_0)+2, ncol = 1,byrow = FALSE,dimnames = list(c(params_txt,params_txt2),"EsT")))
    cat("\n")

    return(list(tableout1,c(params_0,params_1),tableout2,plotout))
  }

}
