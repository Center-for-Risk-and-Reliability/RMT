# Bayesian Step-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2022

stepstress.BAYESest <- function(pt_est,data,stepstresstable,ls,dist,confid,priors,nsamples,burnin,nchains=4){
  #Load pracma library for erf
  library(pracma)
  library(StanHeaders)
  library(rstan)
  library(ggplot2)
  library(shinystan)
  library(cmdstanr)
  library(bayesplot)

  # Add input to this to include prior estimates for LS parameters.
  # Example: priors<-c("normal(3,4),normal(1,4), lognormal(-2,3)")
  # I will have to cite the Rstan text for distributions in the code.  Use lookup("") for the translation.
  # The code takes these and separates them so that they are written into the stan file.
  # Then the code will run the program and compute the Bayes estimation

  # Check to see if dist="Exponential" so you can exclude life
  # distribution parameters.
  if (dist=="Exponential") {
    ishift<-0
  } else {
    ishift<-1
  }
  # Check to see if confidence exists
  if(missing(confid)){
    conf.level <- 0.95
  } else {
    conf.level <- confid
  }

  # Resort the input data and then separate censored from non-censored data
  stpstrdatsort<-stepstress.data.cum(data,stepstresstable)
  xirctot<-sort.xircstressdata(stpstrdatsort[[1]])

  # The failure and survival stress levels are assigned a vector of the same
  # sizes as the adjusted time vectors for the Bayes operation.
  SF <- xirctot[[3]]
  Sc <- xirctot[[4]]
  SFn <- rep(last(stepstresstable[,1]),length(SF))
  Scn <- rep(last(stepstresstable[,1]),length(Sc))

  # State or pull TTF, TTS, tvecti, and tvectj for the tau vector functions
  TTF <- xirctot[[1]]
  TTS <- xirctot[[2]]
  tvecti <- stpstrdatsort[[3]]
  tvectj <- stpstrdatsort[[4]]

  # Obtain the cumulative sum of the end times
  t_end <- cumsum(stepstresstable[,dim(stepstresstable)[[2]]])

  # Compute number of stresses
  Nstress<-dim(stepstresstable)[2]-1

  # Pull the stresses from the table for tau evaluation
  S1i <- stepstresstable[,1]

  if(Nstress==2){
    S2i <- stepstresstable[,2]
  }

  if(Nstress>=3){
    Si <- stepstresstable[,1:Nstress]
  }

  # Setup the delta-t i and j matrices, the single stress i and j matrices, and
  # the single stress i and j vectors
  # Note that this is only for a single stress case and also only for Weibull,
  # Exponential, and Lognormal distributions.  Normal and 2P Exponential and
  # multi-stress cases to come.
  Tendv_i <- rep(0,length(TTF))
  Tendv_j <- rep(0,length(TTS))
  delt_i <- matrix(rep(0,length(stpstrdatsort[[6]])*(max(stpstrdatsort[[6]])-1)),nrow = length(stpstrdatsort[[6]]),ncol = max(stpstrdatsort[[6]])-1,byrow = FALSE)
  delt_j <- matrix(rep(0,length(stpstrdatsort[[7]])*(max(stpstrdatsort[[7]])-1)),nrow = length(stpstrdatsort[[7]]),ncol = max(stpstrdatsort[[7]])-1,byrow = FALSE)
  Spartv_i <- rep(0,length(TTF))
  Spartv_j <- rep(0,length(TTS))
  Spartv_i1 <- rep(stpstrdatsort[[5]][1,1],length(TTF))
  Spartv_j1 <- rep(stpstrdatsort[[5]][1,1],length(TTS))
  Sfull_m_num_i <-delt_i
  Sfull_m_denom_i <-delt_i
  Sfull_m_num_j <-delt_j
  Sfull_m_denom_j <-delt_j
  Svec_num_i<-rep(0,max(stpstrdatsort[[6]])-1)
  Svec_denom_i<-rep(0,max(stpstrdatsort[[6]])-1)
  Svec_num_j<-rep(0,max(stpstrdatsort[[7]])-1)
  Svec_denom_j<-rep(0,max(stpstrdatsort[[7]])-1)

  for(i in 1:max(stpstrdatsort[[6]])-1){
    if(i==1){
      delt_i[which(stpstrdatsort[[6]]>=i+1),i]<-t_end[i]
    } else {
      delt_i[which(stpstrdatsort[[6]]>=i+1),i]<-t_end[i]-t_end[i-1]
    }
    Sfull_m_denom_i[,i]<-stpstrdatsort[[5]][i+1,1]
    Svec_num_i[i]<-stpstrdatsort[[5]][i+1,1]
    Svec_denom_i[i]<-stpstrdatsort[[5]][i,1]
  }
  for(i in 1:max(stpstrdatsort[[6]])){
    Sfull_m_num_i[which(stpstrdatsort[[6]]==i),]<-stpstrdatsort[[5]][i,1]
    Spartv_i[which(stpstrdatsort[[6]]==i)]<-stpstrdatsort[[5]][i,1]
  }
  for(i in 1:max(stpstrdatsort[[7]])-1){
    if(i==1){
      delt_j[which(stpstrdatsort[[7]]>=i+1),i]<-t_end[i]
    } else {
      delt_j[which(stpstrdatsort[[7]]>=i+1),i]<-t_end[i]-t_end[i-1]
    }
    Sfull_m_denom_j[,i]<-stpstrdatsort[[5]][i+1,1]
    Svec_num_j[i]<-stpstrdatsort[[5]][i+1,1]
    Svec_denom_j[i]<-stpstrdatsort[[5]][i,1]
  }
  for(i in 1:max(stpstrdatsort[[7]])){
    Sfull_m_num_j[which(stpstrdatsort[[7]]==i),]<-stpstrdatsort[[5]][i,1]
    Spartv_j[which(stpstrdatsort[[7]]==i)]<-stpstrdatsort[[5]][i,1]
  }

  # Check to see if estimate exists
  if(missing(pt_est)){
    pt_est <- 'random'
  } else {
    pt_est <- pt_est
  }

  # Check to see if burn-in exists
  if(missing(burnin)){
    burnin <- floor(nsamples/2)
  } else {
    burnin <- burnin
  }

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((a*Si[i] + b)/(a*Sni[i] + b));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((a*Sj[j] + b)/(a*Snj[j] + b));"

    # Life functions
    lifeF <- "Lifei[i] = b + Si[i]*a;"
    loglifeF <- "Lifei[i] = log(b + Si[i]*a);"
    lifeC <- "Lifej[j] = b + Sj[j]*a"
    loglifeC <- "Lifej[j] = log(b + Sj[j]*a);"
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp(a*(Si[i] - Sni[i]));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp(a*(Sj[j] - Snj[j]));"


    # Life functions
    lifeF <- "Lifei[i] = b*exp(a*Si[i]);"
    loglifeF <- "Lifei[i] = log(b) + a*Si[i];"
    lifeC <- "Lifej[j] = b*exp(a*Sj[j]);"
    loglifeC <- "Lifej[j] = log(b) + a*Sj[j];"
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    lsparams <- "real<lower=0> Ea; real<lower=0> b;"
    lsparamsvec <- c("Ea","b")
    pr1<-paste(c("Ea ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = b*exp(Ea/((8.617385e-5)*Si[i]));"
    loglifeF <- "Lifei[i] = log(b) + (Ea/((8.617385e-5)*Si[i]));"
    lifeC <- "Lifej[j] = b*exp(Ea/((8.617385e-5)*Sj[j]));"
    loglifeC <- "Lifej[j] = log(b) + (Ea/((8.617385e-5)*Sj[j]));"
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = (b/Si[i])*exp(a/Si[i]);"
    loglifeF <- "Lifei[i] = log(b) - log(Si[i]) + (a/Si[i]);"
    lifeC <- "Lifej[j] = (b/Sj[j])*exp(a/Sj[j]);"
    loglifeC <- "Lifej[j] = log(b) - log(Sj[j]) + (a/Sj[j]);"
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = (1/Si[i])*exp(-(a - (b/Si[i])));"
    loglifeF <- "Lifei[i] = -log(Si[i]) - a + (b/Si[i]);"
    lifeC <- "Lifej[j] = (1/Sj[j])*exp(-(a - (b/Sj[j])));"
    loglifeC <- "Lifej[j] = -log(Sj[j]) - a + (b/Sj[j]);"
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((Si[i]/Sni[i])^a);"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((Sj[j]/Snj[j])^a);"

    # Life functions
    lifeF <- "Lifei[i] = b*(Si[i]^a);"
    loglifeF <- "Lifei[i] = log(b) + a*log(Si[i]);"
    lifeC <- "Lifej[j] = b*(Sj[j]^a);"
    loglifeC <- "Lifej[j] = log(b) + a*log(Sj[j]);"
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real<lower=0> b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((Sni[i]/Si[i])^a);"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((Snj[j]/Sj[j])^a);"

    # Life functions
    lifeF <- "Lifei[i] = b *(Si[i]^-a);"
    loglifeF <- "Lifei[i] = log(b) - a*log(Si[i]);"
    lifeC <- "Lifej[j] = b*(Sj[j]^-a);"
    loglifeC <- "Lifej[j] = log(b) - a*log(Sj[j]);"
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    lsparams <- "real a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*((a*log(Si[i]) + b)/(a*log(Sni[i]) + b));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*((a*log(Sj[j]) + b)/(a*log(Snj[j]) + b));"

    # Life functions
    lifeF <- "Lifei[i] = a*log(Si[i]) + b;"
    loglifeF <- "Lifei[i] = log(a*log(Si[i]) + b);"
    lifeC <- "Lifej[j] = a*log(Sj[j]) + b;"
    loglifeC <- "Lifej[j] = log(a*log(Sj[j]) + b);"
  }

  if (ls=="MultiStress") {
    # CHECK THIS LAST
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    lsparams <- "vector a; real b;"
    lsparamsvec <- c("a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2),collapse = " ")

    lifeF <- exp(theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF))
    function(theta) {
      exp(theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF))
    }

    loglifeF <- function(theta) {
      theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF)
    }
    if(missing(Tc)==FALSE){
      lifeC <- function(theta) {
        exp(theta[ishift+1:length(Sc)+ishift+1]%*%c(1,Sc))
      }
      loglifeC <- function(theta) {
        theta[ishift+1:length(Sc)+ishift+1]%*%c(1,Sc)
      }
    }
  }

  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    lsparams <- "real A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("A ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("a ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("b ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "Lifei[i] = A*exp((a/Sf[,1]) + (b/Sf[,2]));"
    loglifeF <- "Lifei[i] = log(A) + (a/Sf[,1]) + (b/Sf[,2]);"
    lifeC <- "Lifej[j] = A*exp((a/Sc[,1]) + (b/Sc[,2]));"
    loglifeC <- "Lifej[j] = log(A) + (a/Sc[,1]) + (b/Sc[,2]);"
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsparams <- "real A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "c/((Sf[,2]^b)*exp(-a/Sf[,1]));"
    loglifeF <- "log(c) - b*log(Sf[,2]) + (a/Sf[,1]);"
    lifeC <- "c/((Sc[,2]^b)*exp(-a/Sc[,1]));"
    loglifeC <- "log(c) - b*log(Sc[,2]) + (a/Sc[,1]);"
  }

  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    lsparams <- "real a; real b; real c; real d"
    lsparamsvec <- c("a","b","c","d")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    pr4<-paste(c("d ~ ",priors[ishift+4],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3,pr4),collapse = " ")

    # Adjusted Times
    Ti_adj <- "tiadj[i] = (Ti[i] - Tendi[i])*exp((Ea/8.617385e-5)*((1/Si[i]) - (1/Sni[i])));"
    Tj_adj <- "tjadj[j] = (Tj[j] - Tendj[j])*exp((Ea/8.617385e-5)*((1/Sj[j]) - (1/Snj[j])));"

    # Life functions
    lifeF <- "(1/Sf[,1])*exp((a + (b/Sf[,1])) + (c + (d/Sf[,1]))*Sf[,2]);"
    loglifeF <- "-log(Sf[,1]) + a + (b/Sf[,1]) + (c + (d/Sf[,1]))*Sf[,2];"
    lifeC <- "(1/Sc[,1])*exp((a + (b/Sc[,1])) + (c + (d/Sc[,1]))*Sc[,2]);"
    loglifeC <- "-log(Sc[,1]) + a + (b/Sc[,1]) + (c + (d/Sc[,1]))*Sc[,2];"
  }

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    distparam <-"real<lower=0> beta;"
    distpriors<-paste(c("beta ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += weibull_lpdf(tiadj | beta, Lifei) + weibull_lccdf(tjadj | beta, Lifej);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("beta",lsparamsvec)
    outputparamset <- c("\U03B2",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")

  }
  if (dist=="Lognormal") {
    distparam <-"real<lower=0> sigma_t;"
    distpriors<-paste(c("sigma_t ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += lognormal_lpdf(tiadj | Lifei, sigma_t) + lognormal_lccdf(tjadj | Lifej, sigma_t);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma_t",lsparamsvec)
    outputparamset <- c("\U03C3_t",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Normal") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += normal_lpdf(tiadj | Lifei, sigma) + normal_lccdf(tjadj | Lifej, sigma);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    outputparamset <- c("\U03C3",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }
  if (dist=="Exponential") {
    loglik <- paste(c("target += exponential_lpdf(tiadj | 1/Lifei) + exponential_lccdf(tjadj 1/Lifej);"),collapse = "")
    params <- lsparams
    paramsvec <- sparamsvec
    outputparamset <- c(lsparamsvec)
    priors <- lspriors
  }
  if (dist=="2PExponential") {
    distparam <-"real<lower=0> sigma;"
    distpriors<-paste(c("sigma ~ ",priors[ishift],";"),collapse = "")

    loglik <- paste(c("target += double_exponential_lpdf(tiadj | Lifei, sigma) + double_exponential_lccdf(tjadj | Lifej, sigma);"),collapse = "")
    params <- paste(c(distparam,lsparams),collapse = " ")
    paramsvec <- c("sigma",lsparamsvec)
    outputparamset <- c("\U03C3",lsparamsvec)
    priors <- paste(c(distpriors,lspriors),collapse = " ")
  }

  # Define stancode here
  # Data Block
  # block1 <- paste(c("data {int<lower=0> n; int<lower=0> m; vector[n] Ti; vector[m] Tj; vector[n] Tendi; vector[m] Tendj; vector[n] Si; vector[m] Sj; vector[n] Sni; vector[m] Snj; vector[n] tiadj; vector[m] tjadj; vector[n] Lifei; vector[m] Lifej; }"),collapse = " ")
  block1 <- paste(c("data {int<lower=0> n; int<lower=0> m; vector[n] Ti; vector[m] Tj; vector[n] Tendi; vector[m] Tendj; vector[n] Si; vector[m] Sj; vector[n] Sni; vector[m] Snj; }"),collapse = " ")
  # Parameter Block
  block2 <- paste(c("parameters {",params,"}"),collapse = " ")
  # Model Block
  if (dist=="Lognormal"){
    block3 <- paste(c("model {vector[n] tiadj; vector[m] tjadj; vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",Ti_adj,loglifeF,"} for(j in 1:m){",Tj_adj,loglifeC,"}",loglik,"}"),collapse = " ")
  } else{
    block3 <- paste(c("model {vector[n] tiadj; vector[m] tjadj; vector[n] Lifei; vector[m] Lifej; ",priors," for(i in 1:n){",Ti_adj,lifeF,"} for(j in 1:m){",Tj_adj,lifeC,"}",loglik,"}"),collapse = " ")
  }
  # Data Input Block
  #datablock <- list(n = length(TTF), m = length(TTS), Ti = TTF, Tj = TTS, Tendi = tvecti, Tendj = tvectj, Si = SF, Sj = Sc, Sni = SFn, Snj = Scn, tiadj = rep(0,length(TTF)), tjadj = rep(0,length(TTS)), Lifei = rep(0,length(TTF)), Lifej = rep(0,length(TTS)))
  datablock <- list(n = length(TTF), m = length(TTS), Ti = TTF, Tj = TTS, Tendi = tvecti, Tendj = tvectj, Si = SF, Sj = Sc, Sni = SFn, Snj = Scn)
  # NOT RUN {
  stanlscode <- paste(c(block1,block2,block3),collapse=" ")
  stanlsfile <- write_stan_file(stanlscode)
  print(stanlsfile)
  # Generate initial list (one list per chain)
  names(pt_est) <- paramsvec
  pt_estlist <- as.list(pt_est)
  init_pt_est <- vector("list",nchains)
  for(i in 1:nchains){
    init_pt_est[[i]] <- pt_estlist
  }
  # Build or compile Stan code to C++
  # return(list(stanlscode,stanlsfile))

  # lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
  lsmod <- cmdstan_model(stanlsfile)
  # return(lsmod)
  # fit <- sampling(lsmod, data = datablock, iter = nsamples, warmup = burnin, init = pt_est)
  fit <- lsmod$sample(data = datablock, init = init_pt_est, chains = nchains, iter_warmup = burnin, iter_sampling = nsamples)
  # }
  # return(fit)
  # Print results.  I need to get this as an output
  # stats <- print(fit, pars = paramsvec, probs=c((1-confid)/2,.5,1-(1-confid)/2))
  # dataout <- fit@.MISC[["summary"]][["msd"]]

  stats <- fit$summary(variables = paramsvec)
  dataout <- fit$draws(format = "df")
  confidbounds <- mcmc_intervals_data(fit$draws(variables = paramsvec),prob_outer = confid)
  outputtable <- matrix(c(stats[[2]],stats[[4]],confidbounds[[5]],stats[[3]],confidbounds[[9]],stats[[8]]), nrow = length(outputparamset), ncol = 6, byrow = FALSE,dimnames = list(outputparamset,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))


  # Trace the Markov Chains for each parameter
  # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  # plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  # plot2_hist <- stan_hist(fit)
  # plot3_density <- stan_dens(fit)
  plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec))
  plot2_hist <- mcmc_hist(fit$draws(paramsvec))
  plot3_density <- mcmc_dens(fit$draws(paramsvec))
  plot4_densityoverlay <- mcmc_dens_overlay(fit$draws(paramsvec))

  # Produce some output text that summarizes the results
  cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
  print(outputtable)
  cat(c("\n"),sep = "")


  return(list(fit,stats,dataout,plot1_MCtrace,plot2_hist,plot3_density,plot4_densityoverlay))

  # stats <- fit$summary(variables = paramsvec)
  # dataout <- fit$draws(format = "df")
  #
  # # Trace the Markov Chains for each parameter
  # # plot1_MCtrace <- traceplot(fit, pars = paramsvec, inc_warmup = TRUE, nrow = 3)
  # # plot1_MCtrace <- mcmc_trace(as.matrix(fit),pars=paramsvec, facet_args = list(nrow = length(paramsvec), labeller = label_parsed))
  # # plot2_hist <- stan_hist(fit)
  # # plot3_density <- stan_dens(fit)
  # plot1_MCtrace <- mcmc_trace(fit$draws(paramsvec))
  # plot2_hist <- mcmc_hist(fit$draws(paramsvec))
  # plot3_density <- mcmc_dens(fit$draws(paramsvec))
  # plot4_densityoverlay <- mcmc_dens_overlay(fit$draws(paramsvec))
  #
  #
  # return(list(fit,stats,dataout,plot1_MCtrace,plot2_hist,plot3_density,plot4_densityoverlay))

}
