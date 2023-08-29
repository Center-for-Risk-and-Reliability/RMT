# Posterior Life-Stress Relationship Plot
# Developed by Dr. Reuel Smith, 2021-2022

post_lifestress_life_plot <- function(post_params,ls,confid,Xlab,Ylab,Smin,Smax,S0){
  library(pracma)
  library(ggplot2)

  # post_params are currently set to the full fit (output[[1]]) until
  # I figure out what the problem is.  Comment out old material until
  # further research can be done to address it.

  # Check to see if confidence exists
  if(missing(confid)){
    conf.level <- 0.95
  } else {
    conf.level <- confid
  }
  #  Check to see if X and Y labels exists
  if(missing(Xlab)){
    Xlab <- "Stress"
  } else {
    Xlab <- Xlab
  }
  if(missing(Ylab)){
    Ylab <- "Life"
  } else {
    Ylab <- Ylab
  }

  # N<-length(post_params[[1]])
  N<-length(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")])
  M<-100
  # Form matrix for Stress ranges
  S_v<-linspace(Smin,Smax,n=M)
  S_m<-matrix(rep(S_v,N),nrow=N,ncol=M,byrow = TRUE)
  if(missing(S0)==FALSE){
    S0_v <- rep(S0,N)
    # S0_v <- matrix(rep(S0,N),nrow=N,ncol=M,byrow = TRUE)
  }

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m + S_m*a_m
    if(missing(S0)==FALSE){
      # life_S0 <- post_params$b + S0_v*post_params$a
      life_S0 <- post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")] + S0_v*post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")]
    }
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*exp(a_m*S_m)
    if(missing(S0)==FALSE){
      # life_S0 <- post_params$b*exp(post_params$a*S0_v)
      life_S0 <- post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]*exp(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")]*S0_v)
    }
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    # Ea_m <- matrix(rep(post_params$Ea,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    Ea_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="Ea")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*exp(Ea_m/((8.617385e-5)*S_m))
    if(missing(S0)==FALSE){
      # life_S0 <- post_params$b*exp(post_params$Ea/((8.617385e-5)*S0_v))
      life_S0 <- post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]*exp(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="Ea")]/((8.617385e-5)*S0_v))
    }
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- (b_m/S_m)*exp(a_m/S_m)
    if(missing(S0)==FALSE){
      # life_S0 <- (post_params$b/S0_v)*exp(post_params$a/S0_v)
      life_S0 <- (post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]/S0_v)*exp(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")]/S0_v)
    }
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- (1/S_m)*exp(-(a_m - (b_m/S_m)))
    if(missing(S0)==FALSE){
      # life_S0 <- (1/S0_v)*exp(-(post_params$a - (post_params$b/S0_v)))
      life_S0 <- (1/S0_v)*exp(-(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")] - (post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]/S0_v)))
    }
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*(S_m^a_m)
    if(missing(S0)==FALSE){
      # life_S0 <- post_params$b*(S0_v^post_params$a)
      life_S0 <- post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]*(S0_v^post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")])
    }
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*(S_m^-a_m)
    if(missing(S0)==FALSE){
      # life_S0 <- post_params$b*(S0_v^-post_params$a)
      life_S0 <- post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]*(S0_v^-post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")])
    }
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- a_m*log(S_m) + b_m
    if(missing(S0)==FALSE){
      # life_S0 <- post_params$a*log(S0_v) + post_params$b
      life_S0 <- post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")]*log(S0_v) + post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")]

    }
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
    # A_m <- matrix(rep(post_params$A,M),nrow=N,ncol=M,byrow = FALSE)
    # a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    # b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    A_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="A")],M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="a")],M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params[[6]][["data"]][["value"]][which(post_params[[6]][["data"]][["parameter"]]=="b")],M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- A_m*exp((a_m/S_m[,1]) + (b_m/S_m[,2]))
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    lsparams <- "real A; real a; real b;"
    lsparamsvec <- c("A","a","b")
    pr1<-paste(c("a ~ ",priors[ishift+1],";"),collapse = "")
    pr2<-paste(c("b ~ ",priors[ishift+2],";"),collapse = "")
    pr3<-paste(c("c ~ ",priors[ishift+3],";"),collapse = "")
    lspriors <- paste(c(pr1,pr2,pr3),collapse = " ")

    lifeF <- "c/((Sf[,2]^b)*exp(-a/Sf[,1]))"
    loglifeF <- "log(c) - b*log(Sf[,2]) + (a/Sf[,1])"
    if(missing(Tc)==FALSE){
      lifeC <- "c/((Sc[,2]^b)*exp(-a/Sc[,1]))"
      loglifeC <- "log(c) - b*log(Sc[,2]) + (a/Sc[,1])"
    }
    if(missing(S0)==FALSE){
      lifepost <- "b + S0ref*a"
      loglifepost <- "log(b + S0ref*a)"
    }
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

    lifeF <- "(1/Sf[,1])*exp((a + (b/Sf[,1])) + (c + (d/Sf[,1]))*Sf[,2])"
    loglifeF <- "-log(Sf[,1]) + a + (b/Sf[,1]) + (c + (d/Sf[,1]))*Sf[,2]"
    if(missing(Tc)==FALSE){
      lifeC <- "(1/Sc[,1])*exp((a + (b/Sc[,1])) + (c + (d/Sc[,1]))*Sc[,2])"
      loglifeC <- "-log(Sc[,1]) + a + (b/Sc[,1]) + (c + (d/Sc[,1]))*Sc[,2]"
    }
    if(missing(S0)==FALSE){
      lifepost <- "(1/S0ref[,1])*exp((a + (b/S0ref[,1])) + (c + (d/S0ref[,1]))*S0ref[,2])"
      loglifepost <- "-log(S0ref[,1]) + a + (b/S0ref[,1]) + (c + (d/S0ref[,1]))*S0ref[,2]"
    }
  }

  # Compute the mean and confidence intervals for the life curves
  life_mean<-colMeans(life_m)
  life_low<-colQuantiles(life_m,(1-confid)/2)
  life_high<-colQuantiles(life_m,1-(1-confid)/2)

  if(missing(S0)==FALSE){
    df<-data.frame(stress=S_v, lifemean = life_mean, lifemin=life_low,lifemax=life_high,life_at_S0=life_S0)
  } else{
    df<-data.frame(stress=S_v, lifemean = life_mean, lifemin=life_low,lifemax=life_high)
  }

  plotout<-ggplot(data=df, aes(stress,lifemean))+ geom_line(color = "blue", size = 0.9)+ geom_ribbon(aes(ymin=lifemin,ymax=lifemax,x=stress), alpha=0.5,fill = "blue")+xlab(Xlab)+ylab(Ylab)
  if(missing(S0)==FALSE){
    stats_life_S0 <- c(mean(life_S0),std(life_S0),quantile(life_S0,c((1-confid)/2,0.5,1-(1-confid)/2)))
    plotout2 <- ggplot(data=df,aes(life_at_S0,fill = "blue", color = "blue"))+geom_density(alpha = 0.1)
  }
  if(missing(S0)==FALSE){
    return(list(plotout,plotout2,S_v,life_mean,life_low,life_high,stats_life_S0))
  } else{
    return(list(plotout,S_v,life_mean,life_low,life_high))
  }
}
