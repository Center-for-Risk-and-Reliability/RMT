# Posterior Life-Stress Reliability Plot
# Developed by Dr. Reuel Smith, 2021-2022

post_lifestress_reliability_plot <- function(post_params,ls,dist,confid,Xmin,Xmax,Xlab,S0,t0,R0){
  library(pracma)
  library(ggplot2)
  # Check to see if confidence exists
  if(missing(confid)){
    conf.level <- 0.95
  } else {
    conf.level <- confid
  }
  #  Check to see if X and Y labels exists
  if(missing(Xlab)){
    Xlab <- "X"
  } else {
    Xlab <- Xlab
  }

  N<-length(post_params[[1]])
  M<-100
  # Form matrix for Stress ranges
  X_v<-linspace(Xmin,Xmax,n=M)
  X_m<-matrix(rep(X_v,N),nrow=N,ncol=M,byrow = TRUE)

  # Initialize life-stress parameter estimates for theta
  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m + S0*a_m
    life_v <- post_params$b + S0*post_params$a
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*exp(a_m*S0)
    life_v <- post_params$b*exp(post_params$a*S0)
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HAS to be in Kelvin for this to work
    Ea_m <- matrix(rep(post_params$Ea,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*exp(Ea_m/((8.617385e-5)*S0))
    life_v <- post_params$b*exp(post_params$Ea/((8.617385e-5)*S0))
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- (b_m/S0)*exp(a_m/S0)
    life_v <- (post_params$b/S0)*exp(post_params$a/S0)
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- (1/S0)*exp(-(a_m - (b_m/S0)))
    life_v <- (1/S0)*exp(-(post_params$a - (post_params$b/S0)))
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*(S0^a_m)
    life_v <- post_params$b*(S0^post_params$a)
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- b_m*(S0^-a_m)
    life_v <- post_params$b*(S0^-post_params$a)
  }

  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- a_m*log(S0) + b_m
    life_v <- post_params$a*log(S0) + post_params$b
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
    A_m <- matrix(rep(post_params$A,M),nrow=N,ncol=M,byrow = FALSE)
    a_m <- matrix(rep(post_params$a,M),nrow=N,ncol=M,byrow = FALSE)
    b_m <- matrix(rep(post_params$b,M),nrow=N,ncol=M,byrow = FALSE)
    life_m <- A_m*exp((a_m/S0[,1]) + (b_m/S0[,2]))
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

  # Fit to log-likelihood distributions
  if (dist=="Weibull") {
    beta_m <- matrix(rep(post_params$beta,M),nrow=N,ncol=M,byrow = FALSE)
    R_m <- 1 - pweibull(X_m,life_m,beta_m)
    if(missing(t0)==FALSE){
      R_t0 <- 1 - pweibull(t0,life_v,post_params$beta)
    }
    if(missing(R0)==FALSE){
      t_R0 <- life_v*(-log(R0))^(1/post_params$beta)
    }
  }
  if (dist=="Lognormal") {
    sigma_t_m <- matrix(rep(post_params$sigma_t,M),nrow=N,ncol=M,byrow = FALSE)
    R_m <- 1 - plnorm(X_m,log(life_m),sigma_t_m)
    if(missing(t0)==FALSE){
      R_t0 <- 1 - plnorm(t0,log(life_v),post_params$sigma_t)
    }
    if(missing(R0)==FALSE){
      t_R0 <- exp(log(life_v) + post_params$sigma_t*qnorm(1-R0))
    }
  }
  if (dist=="Normal") {
    sigma_m <- matrix(rep(post_params$sigma,M),nrow=N,ncol=M,byrow = FALSE)
    R_m <- 1 - pnorm(X_m,life_m,sigma_m)
    if(missing(t0)==FALSE){
      R_t0 <- 1 - pnorm(t0,life_v,post_params$sigma)
    }
    if(missing(R0)==FALSE){
      t_R0 <- life_v + post_params$sigma*qnorm(1-R0)
    }
  }
  if (dist=="Exponential") {
    R_m <- 1 - pexp(X_m,1/life_m)
    if(missing(t0)==FALSE){
      R_t0 <- 1 - pexp(t0,1/life_v)
    }
    if(missing(R0)==FALSE){
      t_R0 <- -life_v*log(R0)
    }
  }
  if (dist=="2PExponential") {
    sigma_m <- matrix(rep(post_params$sigma,M),nrow=N,ncol=M,byrow = FALSE)
    R_m <- 1 - pexp(X_m-life_m,1/sigma_m)
    if(missing(t0)==FALSE){
      R_t0 <- 1 - pexp(t0-life_v,1/post_params$sigma)
    }
    if(missing(R0)==FALSE){
      t_R0 <- life_v - post_params$sigma*log(R0)
    }
  }

  # Compute the mean and confidence intervals for the life curves
  R_mean<-colMeans(R_m)
  R_low<-colQuantiles(R_m,(1-confid)/2)
  R_high<-colQuantiles(R_m,1-(1-confid)/2)

  if(missing(t0)==TRUE&&missing(R0)==TRUE){
    df<-data.frame(X=X_v, Rmean = R_mean, Rmin=R_low,Rmax=R_high)
  }
  if(missing(t0)==FALSE&&missing(R0)==TRUE){
    df<-data.frame(X=X_v, Rmean = R_mean, Rmin=R_low,Rmax=R_high,reliability_at_t0=R_t0)
  }
  if(missing(t0)==TRUE&&missing(R0)==FALSE){
    df<-data.frame(X=X_v, Rmean = R_mean, Rmin=R_low,Rmax=R_high,life_at_R0=t_R0)
  }
  if(missing(t0)==FALSE&&missing(R0)==FALSE){
    df<-data.frame(X=X_v, Rmean = R_mean, Rmin=R_low,Rmax=R_high,life_at_R0=t_R0,reliability_at_t0=R_t0)
  }
  plotout<-ggplot(data=df, aes(X,Rmean)) + geom_line(color = "red", size = 0.9) + geom_ribbon(aes(ymin=Rmin,ymax=Rmax,x=X), alpha=0.5,fill = "red")+xlab(Xlab)+ylab("Reliability")

  if(missing(t0)==FALSE){
    stats_R_t0 <- c(mean(R_t0),std(R_t0),quantile(R_t0,c((1-confid)/2,0.5,1-(1-confid)/2)))
    plotout2 <- ggplot(data=df,aes(reliability_at_t0,fill = "blue", color = "blue"))+geom_density(alpha = 0.1)
  }
  if(missing(R0)==FALSE){
    stats_t_R0 <- c(mean(t_R0),std(t_R0),quantile(t_R0,c((1-confid)/2,0.5,1-(1-confid)/2)))
    plotout3 <- ggplot(data=df,aes(life_at_R0,fill = "blue", color = "blue"))+geom_density(alpha = 0.1)
  }

  if(missing(t0)==TRUE&&missing(R0)==TRUE){
    return(list(plotout,X_v,R_mean,R_low,R_high))
  }
  if(missing(t0)==FALSE&&missing(R0)==TRUE){
    return(list(plotout,plotout2,X_v,R_mean,R_low,R_high,stats_R_t0))
  }
  if(missing(t0)==TRUE&&missing(R0)==FALSE){
    return(list(plotout,plotout3,X_v,R_mean,R_low,R_high,stats_t_R0))
  }
  if(missing(t0)==FALSE&&missing(R0)==FALSE){
    return(list(plotout,plotout2,plotout3,X_v,R_mean,R_low,R_high,stats_R_t0,stats_t_R0))
  }
}
