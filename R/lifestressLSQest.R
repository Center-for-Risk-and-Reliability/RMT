# Least-Squares Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2022

lifestress.LSQest <- function(ls,dist,pp,therm) {
  #Load pracma library for pseudoinverse
  library(pracma)

  # First check and see that there are multiple stress levels
  if(length(pp)<3) {
    stop('Need more than one stress level to generate estimates')
  }
  # Then check and see if there are single entry data
  if(length(pp)%%3==0){
    singledat<-0 # FALSE Single data does not exist
  } else{
    singledat<-1 # TRUE Single data exists
  }

  if(length(pp[[1]]==2)){
    if(missing(therm)){
      therm<-1
      alttherm<-2
    } else {
      if(therm==1){
        alttherm<-2
      }
      if(therm==2){
        alttherm<-1
      }
    }
  }

  # Setup vectors (for cases with and without single point data)
  if(singledat==0){
    # Sets up existing probability plot curve life and stress vectors
    L<-rep(0,length(pp)/3)
    if (length(pp[[1]])<2){
      S<-rep(0,length(pp)/3)
    } else {
      S<-matrix(rep(0,(length(pp)/3)*length(pp[[1]])),nrow=length(pp)/3,ncol=length(pp[[1]]),byrow = TRUE)
    }
    distparams<-rep(0,length(pp)/3)
  } else if(singledat==1){
    # Sets up probability plot curve and single entry L-S life and stress vectors
    L<-rep(0,(length(pp)-1)/3 + length(tail(pp,n=1)[[1]]))
    if (length(pp[[1]])<2){
      S<-rep(0,(length(pp)-1)/3 + length(tail(pp,n=1)[[1]]))
    } else {
      # NOTE TEST THIS UNDER APPROPRIATE CIRCUMSTANCES
      S<-matrix(rep(0,((length(pp)-1)/3 + length(tail(pp,n=1)[[1]]))*length(pp[[1]])),nrow=(length(pp)-1)/3 + length(tail(pp,n=1)[[1]]),ncol=length(pp[[1]]),byrow = TRUE)
    }
    # Distribution parameter pulls ONLY apply to the probability plots
    distparams<-rep(0,(length(pp)-1)/3)
  }

  # Fill in Stress and Life Vectors
  if(singledat==0){
    for(i2 in 1:(length(pp)/3)){
      # Stress Levels
      if (length(pp[[1]])<2){
        S[i2]<-pp[[i2*3-2]]
      } else {
        for(j in 1:length(pp[[1]])){
          S[i2,j] <- pp[[i2*3-2]][[j]]
        }
      }

      # Life Estimates
      if (dist=="Weibull") {
        L[i2]<-pp[[i2*3-1]][,1]
        distparams[i2]<-pp[[i2*3-1]][2]
      }
      if (dist=="Lognormal") {
        L[i2]<-exp(pp[[i2*3-1]][,1])
        distparams[i2]<-pp[[i2*3-1]][2]
      }
      if (dist=="Normal") {
        L[i2]<-pp[[i2*3-1]][,1]
        distparams[i2]<-pp[[i2*3-1]][2]
      }
      if (dist=="Exponential") {
        L[i2]<-1/pp[[i2*3-1]][,1]
      }
      if (dist=="2PExponential") {
        L[i2]<-pp[[i2*3-1]][,1]+pp[[i2*3-1]][,2]
        distparams[i2]<-pp[[i2*3-1]][2]
      }
    }
  } else if(singledat==1){
    # First Tabulate Probability Plot S and L data
    for(i2 in 1:((length(pp)-1)/3)){
      # Stress Levels
      if (length(pp[[1]])<2){
        S[i2]<-pp[[i2*3-2]]
      } else {
        for(j in 1:length(pp[[1]])){
          S[i2,j] <- pp[[i2*3-2]][[j]]
        }
      }

      # Life Estimates
      if (dist=="Weibull") {
        L[i2]<-pp[[i2*3-1]][,1]
        distparams[i2]<-pp[[i2*3-1]][2]
      }
      if (dist=="Lognormal") {
        L[i2]<-exp(pp[[i2*3-1]][,1])
        distparams[i2]<-pp[[i2*3-1]][2]
      }
      if (dist=="Normal") {
        L[i2]<-pp[[i2*3-1]][,1]
        distparams[i2]<-pp[[i2*3-1]][2]
      }
      if (dist=="Exponential") {
        L[i2]<-1/pp[[i2*3-1]][,1]
      }
      if (dist=="2PExponential") {
        L[i2]<-pp[[i2*3-1]][,1]+pp[[i2*3-1]][,2]
        distparams[i2]<-pp[[i2*3-1]][2]
      }
    }
    # Next tabulate the single point data
    for(i2 in 1:length(tail(pp,n=1)[[1]])){
      # S[i2+(length(pp)-1)/3]<-tail(pp,n=1)[[1]][[i2]][,3]
      # L[i2+(length(pp)-1)/3]<-tail(pp,n=1)[[1]][[i2]][,1]
      S[i2+(length(pp)-1)/3]<-tail(pp,n=1)[[1]][[i2]][[3]]
      L[i2+(length(pp)-1)/3]<-tail(pp,n=1)[[1]][[i2]][[1]]
    }
  }
  # return(list(S,L))

  if (dist=="Weibull") {
    # Writeup for the output text
    dist_txt<-dist
    distparam_txt<-"\U03B2"
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
    life_txt2<-"b*exp(Ea/(K*S))"
    loglife_txt<-"(log(b) + (Ea/(K*S)))"
  }
  if (ls=="Eyring"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    Lvals<-log(L)
    params  <- nls(Lvals ~ log(b) -log(S) + (a/S),start = list(a = 1,b = 3))
    lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
    R2 <- summary(params)$r.squared
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
    R2 <- summary(params)$r.squared
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
    params_txt<-c("a","b")
    # Writeup for the output text
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
    if(length(pp[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),S),nrow=length(pp)/3,ncol=1+length(pp[[1]]),byrow=FALSE)
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
    if(length(pp[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    Lvals<-log(L)
    Svals<-matrix(c(rep(1,length(S[,1])),1/S[,therm],1/S[,alttherm]),nrow=length(pp)/3,ncol=3,byrow=FALSE)
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
    if(length(pp[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    Lvals<-log(L)
    Svals<-matrix(c(1/S[,therm],-log(S[,alttherm]),rep(1,length(S[,1]))),nrow=length(pp)/3,ncol=3,byrow=FALSE)
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
    if(length(pp[[1]])<2) {
      stop('Select a data set with more than one stress type.')
    }
    Lvals<-log(L)+log(S[,therm])
    Svals<-matrix(c(rep(1,length(S[,1])),1/S[,therm],S[,alttherm],S[,alttherm]/S[,therm]),nrow=length(pp)/3,ncol=4,byrow=FALSE)
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


  # Writeup for the output text
  if (dist=="Weibull") {
    dist_txt<-dist
    distparam_txt<-"\U03B2"
    pdf_txt<-c("(\U03B2/",life_txt2,")*(t/",life_txt2,")^(\U03B2-1)*exp(-(t/",life_txt2,")^\U03B2)")
    life_txt<-"63.2% Life - \U03B1"
  }
  if (dist=="Lognormal") {
    dist_txt<-dist
    distparam_txt<-"\U03C3_t"
    pdf_txt<-c("[1/(\U03C3_t t\U221A 2\U03C0)]exp[-0.5*\U03C3_t^(-2)*(ln(t) - ",loglife_txt,")^2]")
    life_txt<-"Median Life - exp(\U03BC_t)"
  }
  if (dist=="Normal") {
    dist_txt<-dist
    distparam_txt<-"\U03C3"
    pdf_txt<-c("[1/(\U03C3 \U221A 2\U03C0)]exp[-0.5*\U03C3^(-2)*(t - ",life_txt2,")^2]")
    life_txt<-"Median Life - \U03BC"
  }
  if (dist=="Exponential") {
    dist_txt<-dist
    pdf_txt<-c("[1/",life_txt2,"]*exp(-t/",life_txt2,")")
    life_txt<-"Median Life 1/\U03BB"
  }
  if (dist=="2PExponential") {
    dist_txt<-"Two-Parameter Exponential"
    distparam_txt<-"\U03C3"
    life_txt<-"Median Life - \U03BC"
  }

  # Group all parameters
  # Check to see if any distribution parameters were tabulated
  if(dist=="Exponential") {
    LSQ<-lsparams
    params_txt<-params_txt
  }
  else {
    LSQ<-c(mean(distparams[which(is.na(distparams) == FALSE)]),lsparams)
    params_txt<-c(distparam_txt,params_txt)
  }

  # Produce some output text that summariZes the results
  cat(c("Least-Squares estimates for the ",ls_txt,"-",dist_txt," Life-Stress model.\n\nf(t,S) = ",pdf_txt,"\n\n"),sep = "")
  print(matrix(c(LSQ), nrow = 1, ncol = length(LSQ), byrow = TRUE,dimnames = list(c("Life-Stress Parameters"),params_txt)))
  cat("\n")
  if(length(pp[[1]])<2){
    print(matrix(c(S,L), nrow = 2, ncol = length(S), byrow = TRUE, dimnames = list(c("Stress",life_txt))))
  } else{
    print(matrix(c(unlist(S),L), nrow = 1+length(pp[[1]]), ncol = length(pp)/3, byrow = TRUE, dimnames = list(c(names(pp[[1]]),life_txt))))
  }
  cat(c("\nCoefficient of Determination R^2 - ",R2))
  cat("\n")

  # Return parameter list
  return(list(S,L,LSQ,R2))
}
