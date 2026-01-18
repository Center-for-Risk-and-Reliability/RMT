# Least-Squares Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2025

lifestress.LSQest <- function(data,ls,dist,pp="Blom",xlabel1="X",
                              Suse=NULL,Llab=NULL,Slab=NULL,Slab2=NULL,param2 = NULL,therm=1,CDFrangesetting = 1,
                              stressunit1 = NULL, stressunit2 = NULL) {
  #Load pracma library for pseudoinverse
  library(pracma)

  # Pull raw data for SSE and R2 of all data
  L.full <- data[,1][which(data[,2]==1)]
  S.full <- data[,3][which(data[,2]==1)]
  if(dim(data)[2]==4){
    S.full.2 <- data[,4][which(data[,2]==1)]
  }

  # Compute probability plotting output first based on input
  # UPDATE (10/9/2025) - ppoutput adjusted to account for changes to probability plotting output
  # UPDATE (11/7/2023) - Now includes the probability plot for the stress levels
  # UPDATE (11/20/2023) - Now includes Gumbel, Logistic, and Log-logistic life distribution options
  if (dist=="Weibull" && is.null(param2)==TRUE) { # Weibull case where both α and β unknown
    output <- probplot.wbl(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Weibull" && is.null(param2)==FALSE) { # Weibull case where α is unknown but β is given
    output <- probplot.wbl(data,pp,xlabel1,setbeta = param2,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="3PWeibull") {
    output <- probplot.wbl3P(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Lognormal") {
    output <- probplot.logn(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Normal") {
    output <- probplot.nor(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Exponential") {
    output <- probplot.exp(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="2PExponential") {
    output <- probplot.exp2P(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Gumbel") {
    output <- probplot.gumb(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Logistic") {
    output <- probplot.logist(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Loglogistic") {
    output <- probplot.loglogist(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="Gamma") {
    output <- probplot.gam(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  if (dist=="3PGamma" || dist=="GeneralizedGamma") {
    output <- probplot.gam3P(data,pp,xlabel1,CDFrangesetting = CDFrangesetting,stressunit1 = stressunit1,stressunit2 = stressunit2)
  }
  ppoutput <- output$output # Main output of stress and parameter estimates
  plotoutput <- output$prob_plot # Secondary output of probability plot

  # return(ppoutput)
  # First check and see that there are multiple stress levels
  # if(length(ppoutput[[1]]$`Stress Level`)<2) {
  #   stop('Need more than one stress level to generate estimates')
  # }
  # Then check and see if there are single entry data
  if(is.null(ppoutput$singledat)==TRUE){
    singledat<-0 # FALSE Single data does not exist
  } else{
    singledat<-1 # TRUE Single data exists
  }

  if(length(ppoutput[[1]]$`Stress Level`)==2){
    if(therm==1){
      alttherm<-2
    }
    if(therm==2){
      alttherm<-1
    }
  }
  # Setup vectors (for cases with and without single point data)
  if(singledat==0){ # Group S, L, and distparams in the case of no single point data.  This is to set up the relation plot.
    # Sets up existing probability plot curve life and stress vectors
    L<-rep(0,length(ppoutput))
    if (length(ppoutput[[1]]$`Stress Level`)<2){
      S<-rep(0,length(ppoutput))
    } else {
      S<-matrix(rep(0,(length(ppoutput))*length(ppoutput[[1]]$`Stress Level`)),nrow=length(ppoutput),ncol=length(ppoutput[[1]]$`Stress Level`),byrow = TRUE)
    }
    distparams<-rep(0,length(ppoutput))
  } else if(singledat==1){
    # Sets up probability plot curve and single entry L-S life and stress vectors
    L<-rep(0,(length(ppoutput[[1]]$`Stress Level`)-1) + length(tail(ppoutput,n=1)[[1]]))
    if (length(ppoutput[[1]]$`Stress Level`)<2){
      S<-rep(0,(length(ppoutput[[1]]$`Stress Level`)-1) + length(tail(ppoutput[[1]]$`Stress Level`,n=1)[[1]]))
    } else {
      # NOTE TEST THIS UNDER APPROPRIATE CIRCUMSTANCES
      S<-matrix(rep(0,((length(ppoutput[[1]]$`Stress Level`)-1) + length(tail(ppoutput[[1]]$`Stress Level`,n=1)[[1]]))*length(ppoutput[[1]]$`Stress Level`)),nrow=(length(ppoutput[[1]]$`Stress Level`)-1) + length(tail(ppoutput[[1]]$`Stress Level`,n=1)[[1]]),ncol=length(ppoutput[[1]]$`Stress Level`),byrow = TRUE)
    }
    # Distribution parameter pulls ONLY apply to the probability plots
    distparams<-rep(0,(length(ppoutput[[1]]$`Stress Level`)-1))
  }
  # return(list(ppoutput,S,L,distparams))
  # Fill in Stress and Life Vectors
  if(singledat==0){
    for(i2 in 1:(length(ppoutput))){
      # Stress Levels
      if (length(ppoutput[[1]]$`Stress Level`)<2){
        S[i2]<-ppoutput[[i2]]$`Stress Level`
      } else {
        for(j in 1:length(ppoutput[[1]]$`Stress Level`)){
          S[i2,j] <- ppoutput[[i2]]$`Stress Level`[j]
        }
      }

      # Life Estimates
      if (dist=="Weibull") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="3PWeibull") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Lognormal") {
        L[i2]<-exp(ppoutput[[i2]]$`Parameter Estimates`[1])
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Normal") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Exponential") {
        L[i2]<-1/ppoutput[[i2]]$`Parameter Estimates`[1]
      }
      if (dist=="2PExponential") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]+ppoutput[[i2]]$`Parameter Estimates`[2]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Gumbel") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Logistic") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Loglogistic") {
        L[i2]<-exp(ppoutput[[i2]]$`Parameter Estimates`[1])
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Gamma") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]*ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="3PGamma" || dist=="GeneralizedGamma") {
        # mean = β x Γ(α + 1/γ)/Γ(α)
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]*(gamma(ppoutput[[i2]]$`Parameter Estimates`[1] + (1/ppoutput[[i2]]$`Parameter Estimates`[3]))/gamma(ppoutput[[i2]]$`Parameter Estimates`[1]))
      }
    }
  } else if(singledat==1){
    # First Tabulate Probability Plot S and L data
    for(i2 in 1:((length(ppoutput)-1))){
      # Stress Levels
      if (length(ppoutput[[1]]$`Stress Level`)<2){
        S[i2]<-ppoutput[[i2]]$`Stress Level`
      } else {
        for(j in 1:length(ppoutput[[1]]$`Stress Level`)){
          S[i2,j] <- ppoutput[[i2]]$`Stress Level`[j]
        }
      }

      # Life Estimates
      if (dist=="Weibull") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="3PWeibull") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Lognormal") {
        L[i2]<-exp(ppoutput[[i2]]$`Parameter Estimates`[1])
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Normal") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Exponential") {
        L[i2]<-1/ppoutput[[i2]]$`Parameter Estimates`[1]
      }
      if (dist=="2PExponential") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]+ppoutput[[i2]]$`Parameter Estimates`[2]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Gumbel") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Logistic") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Loglogistic") {
        L[i2]<-exp(ppoutput[[i2]]$`Parameter Estimates`[1])
        distparams[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="Gamma") {
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[1]*ppoutput[[i2]]$`Parameter Estimates`[2]
      }
      if (dist=="3PGamma" || dist=="GeneralizedGamma") {
        # mean = β x Γ(α + 1/γ)/Γ(α)
        L[i2]<-ppoutput[[i2]]$`Parameter Estimates`[2]*(gamma(ppoutput[[i2]]$`Parameter Estimates`[1] + (1/ppoutput[[i2]]$`Parameter Estimates`[3]))/gamma(ppoutput[[i2]]$`Parameter Estimates`[1]))
      }
    }
    # Next tabulate the single point data
    for(i2 in 1:length(tail(ppoutput,n=1)[[1]])){
      # Enter NA for distparams of single point data
      S[i2+(length(ppoutput)-1)]<-tail(ppoutput,n=1)[[1]][[i2]][[3]]
      L[i2+(length(ppoutput)-1)]<-tail(ppoutput,n=1)[[1]][[i2]][[1]]
      distparams[i2+(length(ppoutput)-1)]<-rep(NA,length(tail(ppoutput,n=1)[[1]][[i2]][[3]]))
    }
  }

  # return(list(ppoutput,S,L,distparams))

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
  if (dist=="Gamma") {
    # Writeup for the output text
    dist_txt<-dist
  }
  if (dist=="3PGamma" || dist=="GeneralizedGamma") {
    # Writeup for the output text
    dist_txt<-"Generalized Gamma"
  }

  # return(list(S,L))
  # LSQ Estimates for Life-Stress Model
  # Executes the LSQ estimates of life-stress model "ls"
  if (ls=="Linear"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(L ~ poly(S, 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])
    SST <- sum((L.full - mean(lsparams[2] + S.full*lsparams[1]))^2)
    SSE <- sum((L.full - (lsparams[2] + S.full*lsparams[1]))^2)
    R2 <- 1 - (SSE/SST)
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
    SST <- sum((log(L.full) - mean(log(lsparams[2]) + S.full*lsparams[1]))^2)
    SSE <- sum((log(L.full) - (log(lsparams[2]) + S.full*lsparams[1]))^2)
    R2 <- 1 - (SSE/SST)
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
    SST <- sum((log(L.full) - mean(log(lsparams[2]) + lsparams[1]/S.full))^2)
    SSE <- sum((log(L.full) - (log(lsparams[2]) + lsparams[1]/S.full))^2)
    R2 <- 1 - (SSE/SST)
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
    # K<-8.617333262e-5
    params  <- lm(log(L) ~ poly(1/S, 1, raw=TRUE))
    lsparams <- c(K*summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    SST <- sum((log(L.full) - mean(log(lsparams[2]) + (lsparams[1]/(K*S.full))))^2)
    SSE <- sum((log(L.full) - (log(lsparams[2]) + (lsparams[1]/(K*S.full))))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("E_a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"b*exp(E_a/(K*S))"
    loglife_txt<-"(log(b) + (E_a/(K*S)))"
  }
  if (ls=="Eyring"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    init.params  <- lm(log(L) ~ poly(1/S, 1, raw=TRUE))
    Lvals<-log(L)
    params  <- nls(Lvals ~ logb - log(S) + (a/S),start = list(a = summary(init.params)$coefficients[2,1],logb = summary(init.params)$coefficients[1,1]))
    # Form Linear Algebra equation A*x = b to solve for parameters
    lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
    lsparams[2] <- exp(lsparams[2])
    SST <- sum((L.full - mean(log(lsparams[2]) - log(S.full) + (lsparams[1]/S.full)))^2)
    SSE <- sum((L.full - (log(lsparams[2]) - log(S.full) + (lsparams[1]/S.full)))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"(b/S)*exp(a/S)"
    loglife_txt<-"(log(b) - log(S) + (a/S))"
  }
  # return(list(L,S,lsparams,summary(params)))

  if (ls=="Eyring2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    Lvals<-log(L)
    params  <- nls(Lvals ~ -log(S) + (b/S) - a,start = list(a = 1,b = 3))
    lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
    SST <- sum((log(L.full) - mean(-log(S.full) - lsparams[1] + (lsparams[2]/S.full)))^2)
    SSE <- sum((log(L.full) - (-log(S.full) - lsparams[1] + (lsparams[2]/S.full)))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-"Eyring (Type-2)"
    life_txt2<-"(1/S)*exp(-(a - (b/S)))"
    loglife_txt<-"(-log(S) - a + (b/S))"
  }
  if (ls=="Power"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(log(S), 1, raw=TRUE))
    lsparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    SST <- sum((log(L.full) - mean(log(lsparams[2]) + lsparams[1]*log(S.full)))^2)
    SSE <- sum((log(L.full) - (log(lsparams[2]) + lsparams[1]*log(S.full)))^2)
    R2 <- 1 - (SSE/SST)
    # Writeup for the output text
    params_txt<-c("a","b")
    ls_txt<-ls
    life_txt2<-"b*(S^a)"
    loglife_txt<-"(ln(b) + aln(S))"
  }
  if (ls=="PowerwithBias"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - c, lsparams[3] - R^2

    Lvals <- L
    PwBSSE <- function(theta){
      SSE <- sum((Lvals - exp(theta[3]) - exp(theta[2])*(S^theta[1]))^2)
      return(SSE)
    }

    params <- nlminb(c(1,1,1),PwBSSE,hessian=TRUE,lower = c(-Inf,-Inf,-Inf),upper = c(Inf,Inf,Inf))$par
    lsparams <- params
    lsparams[2] <- exp(lsparams[2])
    lsparams[3] <- exp(lsparams[3])
    SST <- sum((L.full - mean((lsparams[2] + lsparams[2]*(S.full^lsparams[1]))))^2)
    SSE <- sum((L.full - (lsparams[3] + lsparams[2]*(S.full^lsparams[1])))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b","c")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"c + b*(S^a)"
    loglife_txt<-"log(c + b*(S^a))"
  }
  if (ls=="InversePower"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    params  <- lm(log(L) ~ poly(log(S), 1, raw=TRUE))
    lsparams <- c(-summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
    SST <- sum((log(L.full) - mean(log(lsparams[2]) - lsparams[1]*log(S.full)))^2)
    SSE <- sum((log(L.full) - (log(lsparams[2]) - lsparams[1]*log(S.full)))^2)
    R2 <- 1 - (SSE/SST)
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
    SST <- sum((log(L.full) - mean(-log(lsparams[2]) - lsparams[1]*log(S.full)))^2)
    SSE <- sum((log(L.full) - (-log(lsparams[2]) - lsparams[1]*log(S.full)))^2)
    R2 <- 1 - (SSE/SST)
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
    SST <- sum((L.full - mean(lsparams[2] + log(S.full)*lsparams[1]))^2)
    SSE <- sum((L.full - (lsparams[2] + log(S.full)*lsparams[1]))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b")
    # Writeup for the output text
    ls_txt<-ls
    life_txt2<-"(b + a*ln(S))"
    loglife_txt<-"ln(b + a*ln(S))"
  }
  if (ls=="MultiStress"){
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    defaultStressname <- c("S\U2081","S\U2082","S\U2083","S\U2084","S\U2085","S\U2086","S\U2087","S\U2088","S\U2089","S\U2081\U2080")
    if(length(ppoutput[[1]]$`Stress Level`)<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]$`Stress Level`))==TRUE){
      names(ppoutput[[1]]$`Stress Level`) <- defaultStressname[1:length(ppoutput[[1]]$`Stress Level`)]
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),S),nrow=length(ppoutput),ncol=1+length(ppoutput[[1]]$`Stress Level`),byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    lnLmodel <- Svals%*%lsparams
    SSE <- sum((Lvals - lnLmodel)^2)
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    params_txt<-paste("a_",c(0:length(S[1,])),sep="")
    # Writeup for the output text
    ls_txt<-"Multi-Stress"
    life_txt2<-"exp(a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n)"
    loglife_txt<-"a_0 + a_1*S_1 + a_2*S_2 + ...+ a_n*S_n"
  }
  if (ls=="TempHumidity"){
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    if(length(ppoutput[[1]]$`Stress Level`)<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]$`Stress Level`))==TRUE){
      names(ppoutput[[1]]$`Stress Level`) <- c("Temperature (K)","RH")
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),1/S[,therm],1/S[,alttherm]),nrow=length(ppoutput),ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    params[1]<-exp(params[1])
    lsparams <- c(params)
    SST <- sum((log(L.full) - mean(log(lsparams[1]) + lsparams[2]/S.full + lsparams[3]/S.full.2))^2)
    SSE <- sum((log(L.full) - (log(lsparams[1]) + lsparams[2]/S.full + lsparams[3]/S.full.2))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("A","a","b")
    # Writeup for the output text
    ls_txt<-"Temperature-Humidity"
    life_txt2<-"A exp(a/S + b/H)"
    loglife_txt<-"ln(A) + a/S + b/H"
  }
  if (ls=="TempNonthermal"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    if(length(ppoutput[[1]]$`Stress Level`)<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]$`Stress Level`))==TRUE){
      names(ppoutput[[1]]$`Stress Level`) <- c("Temperature (K)","Nonthermal Stress")
    }
    Lvals<-log(L)
    Svals<-matrix(c(1/S[,therm],-log(S[,alttherm]),rep(1,length(S[,1]))),nrow=length(ppoutput),ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    params[3]<-exp(params[3])
    lsparams <- c(params)
    SST <- sum((log(L.full) - mean(lsparams[1]*(1/S.full) - lsparams[2]*log(S.full.2) + log(lsparams[3])))^2)
    SSE <- sum((log(L.full) - (lsparams[1]*(1/S.full) - lsparams[2]*log(S.full.2) + log(lsparams[3])))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("a","b","c")
    # Writeup for the output text
    ls_txt<-"Temperature-Non-thermal"
    life_txt2<-"c/(U^b * exp(-a/S))"
    loglife_txt<-"a*(1/S) - b*ln(U) + ln(c)"
  }
  if (ls=="Eyring3"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    if(length(ppoutput[[1]]$`Stress Level`)<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]$`Stress Level`))==TRUE){
      names(ppoutput[[1]]$`Stress Level`) <- c("Temperature (K)","Nonthermal Stress")
    }
    Lvals<-log(L)+log(S[,therm])
    Svals<-matrix(c(rep(1,length(S[,1])),1/S[,therm],S[,alttherm],S[,alttherm]/S[,therm]),nrow=length(ppoutput),ncol=4,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    SST <- sum((log(L.full) - mean(-log(S.full) + (lsparams[1] + (lsparams[2]/S.full)) + (lsparams[3] + (lsparams[4]/S.full))*S.full.2))^2)
    SSE <- sum((log(L.full) - (-log(S.full) + (lsparams[1] + (lsparams[2]/S.full)) + (lsparams[3] + (lsparams[4]/S.full))*S.full.2))^2)
    R2 <- 1 - (SSE/SST)
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
    if(length(ppoutput[[1]]$`Stress Level`)<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]$`Stress Level`) <- c("Temperature (K)","Nonthermal Stress")
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),-log(S[,alttherm]),1/S[,therm]),nrow=length(ppoutput),ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    lsparams[1]<-exp(lsparams[1])
    lsparams[3]<-K*lsparams[3]
    SST <- sum((log(L.full) - mean(log(lsparams[1]) + (lsparams[3]/(K*S.full)) - lsparams[2]*log(S.full.2)))^2)
    SSE <- sum((log(L.full) - (log(lsparams[1]) + (lsparams[3]/(K*S.full)) - lsparams[2]*log(S.full.2)))^2)
    R2 <- 1 - (SSE/SST)
    params_txt<-c("A","b","E_a")
    # Writeup for the output text
    ls_txt<-"Eyring (Type 3)"
    life_txt2<-"A exp(E_a/(K*S)) U^-b"
    loglife_txt<-"ln(A) + (E_a/(K*S)) - b ln(U)"
  }
  if (ls=="PH1"){
    # lsparams[1] - parameter beta_0, lsparams[2] - parameter beta_1, lsparams[3] - parameter beta_2
    if(length(ppoutput[[1]]$`Stress Level`)<2) {
      stop('Select a data set with more than one stress type.')
    }
    if(is.null(names(ppoutput[[1]]))==TRUE){
      names(ppoutput[[1]]$`Stress Level`) <- c("Temperature (K)","RH")
    }

    Lvals<-log(L)
    Svals<-matrix(c(rep(-1,length(S[,1])),-1/S[,therm],-1/S[,alttherm]),nrow=length(ppoutput),ncol=3,byrow=FALSE)
    params  <- pinv(Svals)%*%Lvals
    lsparams <- c(params)
    Lmodel <- exp(-lsparams[1])*exp(-lsparams[2]/S[,therm] - lsparams[3]/S[,alttherm])
    lnLmodel <- log(Lmodel)
    R2 <- 1 - sum((Lvals - lnLmodel)^2)/sum((Lvals - mean(Lvals))^2)
    SSE <- sum((Lvals - lnLmodel)^2)
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
  if (dist=="Gamma") {
    dist_txt<-dist
    pdf_txt<-c("[1/",life_txt2,"]*exp(-x/",life_txt2,")")
    life_txt<-"Mean Life 1/\U03BB"
  }
  if (dist=="3PGamma" || dist=="GeneralizedGamma") {
    dist_txt<-dist
    pdf_txt<-c("[1/",life_txt2,"]*exp(-x/",life_txt2,")")
    life_txt<-"Mean Life 1/\U03BB"
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
  if(dist=="2PExponential" || dist=="Weibull" || dist=="Normal" || dist=="Lognormal" || dist=="Gumbel" || dist=="Logistic" || dist=="Loglogistic") {
    LSQ<-c(mean(distparams[which(is.na(distparams) == FALSE)]),lsparams)
    params_txt<-c(distparam_txt,params_txt)
  }
  # NOTE: Uncomment this line if you don't need the relationship plot and just need the output parameters and probability plot
  # return(list(Stress = S,Rated.Life = L,LSQ.point.estimates = LSQ,R2 = R2,SSE = SSE,plotoutput=plotoutput))

  # ====================================================

  # lifestress.relationplot.LSQ(data,ls,dist,params,S=NULL,L=NULL,Smin=NULL,Smax=NULL,Suse=NULL,therm=1,confid=0.95,Llab="Characteristic Life - X",Slab="Characteristic Stress - S") {
  # return(list(S,L))

  if(ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    if(is.null(Suse)==TRUE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,LSQ,S,L,distparams=distparams,stressunit1 = stressunit1)
    }
    if(is.null(Suse)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,LSQ,S,L,distparams=distparams,Suse = Suse,stressunit1 = stressunit1)
    }
    if(is.null(Llab)==FALSE && is.null(Slab)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,LSQ,S,L,distparams=distparams,Suse = Suse,Llab = Llab,Slab = Slab,stressunit1 = stressunit1)
    }
  }
  if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){ # Have working for RAMS
    # return(list(data,ls,dist,LSQ,S,L,distparams,distparams))
    if(is.null(Suse)==TRUE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,LSQ,S,L,distparams=distparams,stressunit1 = stressunit1,stressunit2 = stressunit2)
    }
    if(is.null(Suse)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,LSQ,S,L,distparams=distparams,Suse = Suse,stressunit1 = stressunit1,stressunit2 = stressunit2)
    }
    if(is.null(Llab)==FALSE && is.null(Slab)==FALSE){
      relplotoutput <- lifestress.relationplot.LSQ.2(data,ls,dist,LSQ,S,L,distparams=distparams,Suse = Suse,Llab = Llab,Slab = Slab,stressunit1 = stressunit1,stressunit2 = stressunit2)
    }
  }
  # return(list(S,L,distparams,relplotoutput))
  # return(relplotoutput)
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
  if(length(ppoutput[[1]]$`Stress Level`)<2){
    print(matrix(c(S,L), nrow = 2, ncol = length(S), byrow = TRUE, dimnames = list(c("Stress",life_txt))))
  } else{
    print(matrix(c(unlist(S),L), nrow = 1+length(ppoutput[[1]]$`Stress Level`), ncol = length(ppoutput), byrow = TRUE, dimnames = list(c(names(ppoutput[[1]]$`Stress Level`),life_txt))))
  }
  cat(c("\nSum of Square Error (SSE) - ",SSE))
  cat(c("\nCoefficient of Determination (R^2) - ",R2))
  cat("\n")

  # Return parameter list
  # return(list(S,L,LSQ,R2,plotoutput=plotoutput))
  # FOR USE WITH RELATION PLOT UPDATE
  if(ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    if(is.null(Suse)==TRUE){
      return(list(Stress = S,Rated.Life = L,LSQ.point.estimates = LSQ,R2 = R2,SSE = SSE,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
    }
    if(is.null(Suse)==FALSE){
      return(list(Stress = S,Rated.Life = L,LSQ.point.estimates = LSQ,R2 = R2,SSE = SSE,Use_Life = relplotoutput$Luse,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
    }
  }
  if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){ # Have working for RAMS
    if(is.null(Suse)==TRUE){
      return(list(Stress = S,Rated.Life = L,LSQ.point.estimates = LSQ, R2 = R2,SSE = SSE,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
    }
    if(is.null(Suse)==FALSE){
      return(list(Stress = S,Rated.Life = L,LSQ.point.estimates = LSQ, R2 = R2,SSE = SSE,Use_Life = relplotoutput$Luse,plotoutput=plotoutput,relplotoutput=relplotoutput$relationplot))
    }
  }
}
