# Least-Squares Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2023

lifestress.LSQest <- function(data,ls,dist,pp,xlabel1="X",therm=1,Suse=NULL,Llab=NULL,Slab=NULL) {
  #Load pracma library for pseudoinverse
  library(pracma)

  # Compute probability plotting output first based on input
  # UPDATE (11/7/2023) - Now includes the probability plot for the stress levels
  # UPDATE (11/20/2023) - Now includes Gumbel, Logistic, and Log-logistic life distribution options
  if (dist=="Weibull") {
    ppoutput <- probplot.wbl(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.wbl(data,pp,xlabel1)$prob_plot
  }
  if (dist=="3PWeibull") {
    ppoutput <- probplot.wbl3P(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.wbl3P(data,pp,xlabel1)$prob_plot
    nonparamoutput <- probplot.wbl3P(data,pp,xlabel1)$summary.nonparametric
  }
  if (dist=="Lognormal") {
    ppoutput <- probplot.logn(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.logn(data,pp,xlabel1)$prob_plot
  }
  if (dist=="Normal") {
    ppoutput <- probplot.nor(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.nor(data,pp,xlabel1)$prob_plot
  }
  if (dist=="Exponential") {
    ppoutput <- probplot.exp(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.exp(data,pp,xlabel1)$prob_plot
  }
  if (dist=="2PExponential") {
    ppoutput <- probplot.exp2P(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.exp2P(data,pp,xlabel1)$prob_plot
  }
  if (dist=="Gumbel") {
    ppoutput <- probplot.gumb(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.gumb(data,pp,xlabel1)$prob_plot
  }
  if (dist=="Logistic") {
    ppoutput <- probplot.logist(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.logist(data,pp,xlabel1)$prob_plot
  }
  if (dist=="Loglogistic") {
    ppoutput <- probplot.loglogist(data,pp,xlabel1)[[1]]
    plotoutput <- probplot.loglogist(data,pp,xlabel1)$prob_plot
  }

  # First check and see that there are multiple stress levels
  if(length(ppoutput)<3) {
    stop('Need more than one stress level to generate estimates')
  }
  # Then check and see if there are single entry data
  if(length(ppoutput)%%3==0){
    singledat<-0 # FALSE Single data does not exist
  } else{
    singledat<-1 # TRUE Single data exists
  }

  if(length(ppoutput[[1]]==2)){
    if(therm==1){
      alttherm<-2
    }
    if(therm==2){
      alttherm<-1
    }
  }

  # Setup vectors (for cases with and without single point data)
  if(singledat==0){
    # Sets up existing probability plot curve life and stress vectors
    L<-rep(0,length(ppoutput)/3)
    if (length(ppoutput[[1]])<2){
      S<-rep(0,length(ppoutput)/3)
    } else {
      S<-matrix(rep(0,(length(ppoutput)/3)*length(ppoutput[[1]])),nrow=length(ppoutput)/3,ncol=length(ppoutput[[1]]),byrow = TRUE)
    }
    distparams<-rep(0,length(ppoutput)/3)
  } else if(singledat==1){
    # Sets up probability plot curve and single entry L-S life and stress vectors
    L<-rep(0,(length(ppoutput)-1)/3 + length(tail(ppoutput,n=1)[[1]]))
    if (length(ppoutput[[1]])<2){
      S<-rep(0,(length(ppoutput)-1)/3 + length(tail(ppoutput,n=1)[[1]]))
    } else {
      # NOTE TEST THIS UNDER APPROPRIATE CIRCUMSTANCES
      S<-matrix(rep(0,((length(ppoutput)-1)/3 + length(tail(ppoutput,n=1)[[1]]))*length(ppoutput[[1]])),nrow=(length(ppoutput)-1)/3 + length(tail(ppoutput,n=1)[[1]]),ncol=length(ppoutput[[1]]),byrow = TRUE)
    }
    # Distribution parameter pulls ONLY apply to the probability plots
    distparams<-rep(0,(length(ppoutput)-1)/3)
  }

  # Fill in Stress and Life Vectors
  if(singledat==0){
    for(i2 in 1:(length(ppoutput)/3)){
      # Stress Levels
      if (length(ppoutput[[1]])<2){
        S[i2]<-ppoutput[[i2*3-2]]
      } else {
        for(j in 1:length(ppoutput[[1]])){
          S[i2,j] <- ppoutput[[i2*3-2]][[j]]
        }
      }

      # Life Estimates
      if (dist=="Weibull") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="3PWeibull") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Lognormal") {
        L[i2]<-exp(ppoutput[[i2*3-1]][,1])
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Normal") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Exponential") {
        L[i2]<-1/ppoutput[[i2*3-1]][,1]
      }
      if (dist=="2PExponential") {
        L[i2]<-ppoutput[[i2*3-1]][,1]+ppoutput[[i2*3-1]][,2]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Gumbel") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Logistic") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Loglogistic") {
        L[i2]<-exp(ppoutput[[i2*3-1]][,1])
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
    }
  } else if(singledat==1){
    # First Tabulate Probability Plot S and L data
    for(i2 in 1:((length(ppoutput)-1)/3)){
      # Stress Levels
      if (length(ppoutput[[1]])<2){
        S[i2]<-ppoutput[[i2*3-2]]
      } else {
        for(j in 1:length(ppoutput[[1]])){
          S[i2,j] <- ppoutput[[i2*3-2]][[j]]
        }
      }

      # Life Estimates
      if (dist=="Weibull") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="3PWeibull") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Lognormal") {
        L[i2]<-exp(ppoutput[[i2*3-1]][,1])
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Normal") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Exponential") {
        L[i2]<-1/ppoutput[[i2*3-1]][,1]
      }
      if (dist=="2PExponential") {
        L[i2]<-ppoutput[[i2*3-1]][,1]+ppoutput[[i2*3-1]][,2]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Gumbel") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Logistic") {
        L[i2]<-ppoutput[[i2*3-1]][,1]
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
      if (dist=="Loglogistic") {
        L[i2]<-exp(ppoutput[[i2*3-1]][,1])
        distparams[i2]<-ppoutput[[i2*3-1]][2]
      }
    }
    # Next tabulate the single point data
    for(i2 in 1:length(tail(ppoutput,n=1)[[1]])){
      # S[i2+(length(ppoutput)-1)/3]<-tail(ppoutput,n=1)[[1]][[i2]][,3]
      # L[i2+(length(ppoutput)-1)/3]<-tail(ppoutput,n=1)[[1]][[i2]][,1]
      S[i2+(length(ppoutput)-1)/3]<-tail(ppoutput,n=1)[[1]][[i2]][[3]]
      L[i2+(length(ppoutput)-1)/3]<-tail(ppoutput,n=1)[[1]][[i2]][[1]]
    }
  }

  # return(list(S,L,distparams))

  if (dist=="Weibull") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03B2"
  }
  if (dist=="3PWeibull") {
    # Writeup for the output text
    dist_txt<-"Three-Parameter Weibull"
    distparam_txt<-c("\U03B2","\U03B3")
  }
  if (dist=="Lognormal") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03C3_t"
  }
  if (dist=="Normal") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03C3"
  }
  if (dist=="Exponential") {
    # Writeup for the output text
    dist_txt<-dist
  }
  if (dist=="2PExponential") {
    # Writeup for the output text
    dist_txt<-"Two-Parameter Exponential"
    distparam_txt<-"\U03C3"
  }
  if (dist=="Gumbel") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03C3"
  }
  if (dist=="Logistic") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03C3"
  }
  if (dist=="Loglogistic") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03C3"
  }

  # LSQ Estimates for Life-Stress Model
  # Executes the LSQ estimates of life-stress model "ls"
  if (ls=="Linear"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(L ~ poly(S, 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])
    R2 <- summary(params)$r.squared
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"(b + S*a)"
    loglife_txt<-"ln(b + S*a)"
  }

  if (ls=="Exponential"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(S, 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    R2 <- summary(params)$r.squared
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"b*exp(a*S)"
    loglife_txt<-"(log(b) + a*S)"
  }
  if (ls=="Exponential2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(1/S, 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    R2 <- summary(params)$r.squared
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"b*exp(a/S)"
    loglife_txt<-"(log(b) + a/S)"
  }
  if (ls=="Arrhenius"){
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b, lsparams[3] - R^2
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    params  <- lm(log(L) ~ poly(1/S, 1, raw=TRUE))
    lsparams <- c(K*summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    R2 <- summary(params)$r.squared
    params_txt<-c("E_a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"b*exp(E_a/(K*S))"
    loglife_txt<-"(log(b) + (E_a/(K*S)))"
  }
  if (ls=="Eyring"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    Lvals<-log(L)
    params  <- nls(Lvals ~ log(b) -log(S) + (a/S),start = list(a = 1,b = 3))
    lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
    SST <- sum((Lvals - mean(Lvals))^2)
    SSE <- deviance(params)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"(b/S)*exp(a/S)"
    loglife_txt<-"(log(b) - log(S) + (a/S))"
  }
  if (ls=="Eyring2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    Lvals<-log(L)
    params  <- nls(Lvals ~ -log(S) + (b/S) - a,start = list(a = 1,b = 3))
    lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
    SST <- sum((Lvals - mean(Lvals))^2)
    SSE <- deviance(params)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-"Eyring (Type-2)"
    life_txt2<-"(1/S)*exp(-(a - (b/S)))"
    loglife_txt<-"(-log(Sf) - a + (b/S))"
  }
  if (ls=="Power"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(log(S), 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    R2 <- summary(params)$r.squared
    # Writeup for the output text
    params_txt<-c("a","b")
    ls_txt<-ls
    life_txt2<-"b*(S^a)"
    loglife_txt<-"(ln(b) + aln(S))"
  }
  if (ls=="InversePower"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(log(S), 1, raw=TRUE))
    lsparams <- c(-summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    R2 <- summary(params)$r.squared
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-"Inverse Power"
    life_txt2<-"b*(S^-a)"
    loglife_txt<-"(ln(b) - aln(S))"
  }
  # UPDATE 11/8/2023: Adding an IPL form of l = 1/(bxS^a)
  if (ls=="InversePower2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(log(S), 1, raw=TRUE))
    lsparams <- c(-summary(params)$coefficients[2,1],exp(-summary(params)$coefficients[1,1]))
    R2 <- summary(params)$r.squared
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-"Inverse Power"
    life_txt2<-"1/[b*(S^a)]"
    loglife_txt<-"(-ln(b) - aln(S))"
  }
  if (ls=="Logarithmic"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(L ~ poly(log(S), 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])
    R2 <- summary(params)$r.squared
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"(b + a*ln(S))"
    loglife_txt<-"ln(b + a*ln(S))"
  }
  if (ls=="MultiStress"){
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    defaultStressname <- c("S\U2081","S\U2082","S\U2083","S\U2084","S\U2085","S\U2086","S\U2087","S\U2088","S\U2089","S\U2081\U2080")
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]) <- defaultStressname[1:length(ppoutput[[1]])]
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),S),nrow=length(ppoutput)/3,ncol=1+length(ppoutput[[1]]),byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    lnLmodel <- Svals%*%lsparams
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    params_txt<-paste("a_",c(0:length(S[1,])),sep="")
    # Writeup for the output text
    ls_txt<-"Multi-Stress"
    life_txt2<-"exp(a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n)"
    loglife_txt<-"a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n"
  }
  if (ls=="TempHumidity"){
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]) <- c("Temperature (K)","RH")
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),1/S[,therm],1/S[,alttherm]),nrow=length(ppoutput)/3,ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    params[1]<-exp(params[1])
    lsparams <- c(params)
    Lmodel <- lsparams[1]*exp(lsparams[2]/S[,therm] + lsparams[3]/S[,alttherm])
    lnLmodel <- log(Lmodel)
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    params_txt<-c("A","a","b")
    # Writeup for the output text
    ls_txt<-"Temperature-Humidity"
    life_txt2<-"A exp(a/S + b/H)"
    loglife_txt<-"ln(A) + a/S + b/H"
  }
  if (ls=="TempNonthermal"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]) <- c("Temperature (K)","Nonthermal Stress")
    }
    Lvals<-log(L)
    Svals<-matrix(c(1/S[,therm],-log(S[,alttherm]),rep(1,length(S[,1]))),nrow=length(ppoutput)/3,ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    params[3]<-exp(params[3])
    lsparams <- c(params)
    lnLmodel <- lsparams[1]*(1/S[,therm]) - lsparams[2]*log(S[,alttherm]) + log(lsparams[3])
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    params_txt<-c("a","b","c")
    # Writeup for the output text
    ls_txt<-"Temperature-Non-thermal"
    life_txt2<-"c/(U^b * exp(-a/S))"
    loglife_txt<-"a(1/S) - b*ln(U) + ln(c)"
  }
  if (ls=="Eyring3"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]) <- c("Temperature (K)","Nonthermal Stress")
    }
    Lvals<-log(L)+log(S[,therm])
    Svals<-matrix(c(rep(1,length(S[,1])),1/S[,therm],S[,alttherm],S[,alttherm]/S[,therm]),nrow=length(ppoutput)/3,ncol=4,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    lnLmodel <- -log(S[,therm]) + lsparams[1] + lsparams[2]/S[,therm] + lsparams[3]*S[,alttherm] + lsparams[4]*(S[,alttherm]/S[,therm])
    R2 <- 1 - sum((log(L) - lnLmodel)^2)/sum((log(L) - mean(log(L)))^2)
    params_txt<-c("a","b","c","d")
    # Writeup for the output text
    ls_txt<-"Eyring (Type 3)"
    life_txt2<-"(1/S) exp((a + (b/S)) + (c + (d/S)) U)"
    loglife_txt<-"-ln(S) + (a + (b/S)) + (c + (d/S)) U"
  }
  if (ls=="Eyring4"){
    # lsparams[1] - parameter A, lsparams[2] - parameter b
    # lsparams[3] - parameter Ea
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]) <- c("Temperature (K)","Nonthermal Stress")
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),-log(S[,alttherm]),1/S[,therm]),nrow=length(ppoutput)/3,ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    lsparams[1]<-exp(lsparams[1])
    lsparams[3]<-K*lsparams[3]
    lnLmodel <- log(lsparams[1]) - lsparams[2]*log(S[,alttherm]) + (lsparams[3]/K)*(1/S[,therm])
    R2 <- 1 - sum((log(L) - lnLmodel)^2)/sum((log(L) - mean(log(L)))^2)
    params_txt<-c("A","b","E_a")
    # Writeup for the output text
    ls_txt<-"Eyring (Type 3)"
    life_txt2<-"A exp(E_a/(K*S)) U^-b"
    loglife_txt<-"ln(A) + (E_a/(K*S)) - b ln(U)"
  }
  if (ls=="PH1"){
    # lsparams[1] - parameter beta_0, lsparams[2] - parameter beta_1, lsparams[3] - parameter beta_2
    if(length(ppoutput[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]) <- c("Temperature (K)","RH")
    }

    Lvals<-log(L)
    Svals<-matrix(c(rep(-1,length(S[,1])),-1/S[,therm],-1/S[,alttherm]),nrow=length(ppoutput)/3,ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    Lmodel <- exp(-lsparams[1])*exp(-lsparams[2]/S[,therm] - lsparams[3]/S[,alttherm])
    lnLmodel <- log(Lmodel)
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    params_txt<-c("\U03B2\U2080","\U03B2\U2081","\U03B2\U2082")
    # Writeup for the output text
    ls_txt<-"Proportional-Hazard"
    life_txt2<-"exp(-\U03B2\U2080) exp(-\U03B2\U2081/S - \U03B2\U2082/RH)"
    loglife_txt<-"-\U03B2\U2080 - \U03B2\U2081/S - \U03B2\U2082/RH"
  }

  # Writeup for the output text
  if (dist=="Weibull") {
    dist_txt<-dist
    distparam_txt<-"\U03B2"
    pdf_txt<-c("(\U03B2/",life_txt2,")*(x/",life_txt2,")^(\U03B2-1)*exp(-(x/",life_txt2,")^\U03B2)")
    life_txt<-"63.2% Life - \U03B1"
  }
  if (dist=="3PWeibull") {
    dist_txt<-"Three-Parameter Weibull"
    distparam_txt<-c("\U03B2","\U03B3")
    pdf_txt<-c("(\U03B2/",life_txt2,")*((x-\U03B3)/",life_txt2,")^(\U03B2-1)*exp(-((x-\U03B3)/",life_txt2,")^\U03B2)")
    life_txt<-"63.2% Life - \U03B1"
  }
  if (dist=="Lognormal") {
    dist_txt<-dist
    distparam_txt<-"\U03C3_t"
    pdf_txt<-c("[1/(\U03C3_t x\U221A 2\U03C0)]exp[-0.5*\U03C3_t^(-2)*(ln(x) - ",loglife_txt,")^2]")
    life_txt<-"Median Life - exp(\U03BC_x)"
  }
  if (dist=="Normal") {
    dist_txt<-dist
    distparam_txt<-"\U03C3"
    pdf_txt<-c("[1/(\U03C3 \U221A 2\U03C0)]exp[-0.5*\U03C3^(-2)*(x - ",life_txt2,")^2]")
    life_txt<-"Mean/Median Life - \U03BC"
  }
  if (dist=="Exponential") {
    dist_txt<-dist
    pdf_txt<-c("[1/",life_txt2,"]*exp(-x/",life_txt2,")")
    life_txt<-"Mean Life 1/\U03BB"
  }
  if (dist=="2PExponential") {
    dist_txt<-"Two-Parameter Exponential"
    distparam_txt<-"\U03C3"
    life_txt<-"Median Life - \U03BC"
  }
  if (dist=="Gumbel") {
    dist_txt<-dist
    distparam_txt<-"\U03C3"
    pdf_txt<-c("(1/\U03C3) exp{(x - ",life_txt2,")/\U03C3 - exp[(x - ",life_txt2,")/\U03C3]}")
    life_txt<-"Mode Life - \U03BC"
  }
  if (dist=="Logistic") {
    dist_txt<-dist
    distparam_txt<-"\U03C3"
    pdf_txt<-c("[1/(\U03C3 \U221A 2\U03C0)]exp[-0.5*\U03C3^(-2)*(x - ",life_txt2,")^2]")
    life_txt<-"Mean/Median Life - \U03BC"
  }
  if (dist=="Loglogistic") {
    dist_txt<-dist
    distparam_txt<-"\U03C3_t"
    pdf_txt<-c("[1/(\U03C3_t x\U221A 2\U03C0)]exp[-0.5*\U03C3_t^(-2)*(ln(x) - ",loglife_txt,")^2]")
    life_txt<-"Median Life - exp(\U03BC_x)"
  }

  # Group all parameters
  # Check to see if any distribution parameters were tabulated
  if(dist=="Exponential") {
    LSQ<-lsparams
    params_txt<-params_txt
  }
  if(dist=="3PWeibull") {
    beta_LSQ<-mean(distparams[which(is.na(distparams) == FALSE)])
    for(i in 1:length(nonparamoutput)){
      if(i == 1){
        xfull <- nonparamoutput[[i]][,1]
        Rfull <- nonparamoutput[[i]][,3]
        alpfull <- rep(L[i],length(nonparamoutput[[i]][,1]))
      } else{
        xfull <- c(xfull,nonparamoutput[[i]][,1])
        Rfull <- c(Rfull,nonparamoutput[[i]][,3])
        alpfull <- c(alpfull,rep(L[i],length(nonparamoutput[[i]][,1])))
      }

    }
    gam_LSQest <- function(gam){
      Fx <- sum((log(-log(Rfull)) - beta_LSQ*log(xfull - gam) + beta_LSQ*log(alpfull))^2)
      return(Fx)
    }
    gamma_LSQ <- nlminb(0.5*min(xfull),gam_LSQest,hessian=TRUE,lower = -Inf,upper = 0.99*min(xfull))$par

    LSQ<-c(beta_LSQ,gamma_LSQ,lsparams)
    params_txt<-c(distparam_txt,params_txt)
  }
  if(dist=="2PExponential" || dist=="Weibull" || dist=="Normal" || dist=="Lognormal" || dist=="Gumbell" || dist=="Logistic" || dist=="Loglogistic") {
    LSQ<-c(mean(distparams[which(is.na(distparams) == FALSE)]),lsparams)
    params_txt<-c(distparam_txt,params_txt)
  }
  # return(list(S,L,LSQ,R2,plotoutput=plotoutput))
  # lifestress.relationplot.LSQ(data,ls,dist,params,S=NULL,L=NULL,Smin=NULL,Smax=NULL,Suse=NULL,therm=1,confid=0.95,Llab="Characteristic Life - X",Slab="Characteristic Stress - S") {
  # return(LSQ)

  if(ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    if(is.null(Suse)==TRUE){
      relplotoutput <- lifestress.relationplot.LSQ(data,ls,dist,LSQ,S,L)
    }
    if(is.null(Suse)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ(data,ls,dist,LSQ,S,L,Suse = Suse)
    }
    if(is.null(Llab)==FALSE && is.null(Slab)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ(data,ls,dist,LSQ,S,L,Suse = Suse,Llab = Llab,Slab = Slab)
    }
  }
  #
  # if(is.null(Suse)==TRUE){
  #   relplotoutput <- lifestress.relationplot.LSQ(data,ls,dist,lsparams,S,L)
  # }
  # if(is.null(Suse)==FALSE){
  #   relplotoutput <- lifestress.relationplot.LSQ(data,ls,dist,lsparams,S,L,Suse = Suse)
  # }

  # return(relplotoutput)

  # return(list(c(names(ppoutput[[1]]),life_txt)))

  # Produce some output text that summarizes the results
  cat(c("Least-Squares estimates for the ",ls_txt,"-",dist_txt," Life-Stress model.\n\nf(x,S) = ",pdf_txt,"\n\n"),sep = "")
  print(matrix(c(LSQ), nrow = 1, ncol = length(LSQ), byrow = TRUE,dimnames = list(c("Life-Stress Parameters"),params_txt)))
  cat("\n")
  if(length(ppoutput[[1]])<2){
    print(matrix(c(S,L), nrow = 2, ncol = length(S), byrow = TRUE, dimnames = list(c("Stress",life_txt))))
  } else{
    print(matrix(c(unlist(S),L), nrow = 1+length(ppoutput[[1]]), ncol = length(ppoutput)/3, byrow = TRUE, dimnames = list(c(names(ppoutput[[1]]),life_txt))))
  }
  cat(c("\nCoefficient of Determination R^2 - ",R2))
  cat("\n")

  # Return parameter list
  # return(list(S,L,LSQ,R2,plotoutput=plotoutput))
  # FOR USE WITH RELATION PLOT UPDATE
  if(ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    if(is.null(Suse)==TRUE){
      return(list(S,L,LSQ,R2,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
    }
    if(is.null(Suse)==FALSE){
      return(list(S,L,LSQ,R2,Use_Life = relplotoutput$Luse,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
    }
  } else{
    return(list(S,L,LSQ,R2,plotoutput=plotoutput))
  }
  #
  # if(is.null(Suse)==TRUE){
  #   return(list(S,L,LSQ,R2,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
  # }
  # if(is.null(Suse)==FALSE){
  #   return(list(S,L,LSQ,R2,Use_Life = relplotoutput$Luse,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
  # }

}
