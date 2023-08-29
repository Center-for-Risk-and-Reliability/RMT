# Least-Squares Step-Stress Estimator
# Developed by Dr. Reuel Smith, 2021-2023

stepstress.LSQest <- function(data,stepstresstable,ls,dist,pp,therm=1) {
  #Load pracma library for pseudo-inverse
  library(pracma)
  library(dplyr)

  # Set which stresses are thermal ans which are non-thermal if you have a dual
  # thermal-nonthermal life-stress model.
  if(therm==1){
    alttherm<-2
  }
  if(therm==2){
    alttherm<-1
  }

  # First sort the data for initial probability plotting after checking to see
  # if it is a single stepstress data/table pair or a set of them
  if(is.list(data[[1]])==FALSE){
    # FOR ONE STEP-STRESS TEST
    stpstrdatsort<-stepstress.data(data,stepstresstable)
    full_stpstrdata <- stpstrdatsort[[1]]
    # Compute the cumulative damage
    cumdmg<-stepstress.damage(data,stepstresstable)
    # Obtain the number of stresses and pertinent steps in this
    Nstress<-dim(stepstresstable)[2]-1
  } else{
    # FOR SET OF STEP-STRESS TESTS TO BE COMPUTED TOGETHER
    stpstrdatsort<-vector(mode = "list", length(data))
    cumdmg<-vector(mode = "list", length(data))
    for(i in 1:length(data)){
      stpstrdatsort[[i]]<-stepstress.data(data[[i]],stepstresstable[[i]])
      if(i == 1){
        full_stpstrdata <- stpstrdatsort[[i]][[1]]
      } else{
        full_stpstrdata <- merge(full_stpstrdata,stpstrdatsort[[i]][[1]], all = TRUE, sort=FALSE)
      }
      # Compute the cumulative damage on test-by-test basis
      cumdmg[[i]]<-stepstress.damage(data[[i]],stepstresstable[[i]])
    }
    colnames(full_stpstrdata) <- colnames(data[[1]], do.NULL = FALSE, prefix = "Obs.")
    # Obtain the number of stresses and pertinent steps in this
    Nstress<-dim(stepstresstable[[1]])[2]-1
  }

  # Tabulate initial LSQ Estimates starting with the distributions
  if (dist=="Weibull") {
    distlifeest <- function(fulldata){
      distoutput <- probplot.wbl(fulldata,pp,colnames(data)[1])[[1]]
      iparam <- -1 + cumsum(rep(1+Nstress+2, ((length(distoutput)-rem(length(distoutput),3))/3)))
      ilife <- -2 + cumsum(rep(1+Nstress+2, ((length(distoutput)-rem(length(distoutput),3))/3)))
      istress <- -3 + cumsum(rep(1+Nstress+2, ((length(distoutput)-rem(length(distoutput),3))/3)))
      lifeest <- unlist(distoutput)[ilife]

      return(list(distoutput,lifeest,istress,ilife,iparam))
    }
  }
  if (dist=="Lognormal") {
    distlifeest <- function(fulldata){
      distoutput <-probplot.logn(fulldata,pp,colnames(data)[1])[[1]]
      lifeest<-exp(unlist(distoutput)[ilife])

      return(list(distoutput,lifeest))
    }
  }
  if (dist=="Normal") {
    distlifeest <- function(fulldata){
      distoutput <-probplot.nor(fulldata,pp,colnames(data)[1])[[1]]
      lifeest<-unlist(distoutput)[ilife]
      return(list(distoutput,lifeest))
    }
  }
  if (dist=="Exponential") {
    distlifeest <- function(fulldata){
      distoutput<-probplot.exp(fulldata,pp,colnames(data)[1])[[1]]
      ilife2 <- c(3,7,11,15)
      # lifeest<-1/unlist(distoutput)[ilife]
      lifeest<-1/unlist(distoutput)[ilife2]
      return(list(distoutput,lifeest))
    }
  }
  if (dist=="2PExponential") {
    distlifeest <- function(fulldata){
      distoutput <-probplot.exp2P(fulldata,pp,colnames(data)[1])[[1]]
      lifeest<-unlist(distoutput)[ilife]
      return(list(distoutput,lifeest))
    }
  }
  output1<-distlifeest(full_stpstrdata)
  distoutput<-output1[[1]]
  lifeest<-output1[[2]]
  # Collect the stresses for evaluation
  if(Nstress==1){
    Stressset<-unlist(distoutput)[output1[[3]]]
  }
  if(Nstress==2){
    Stressset<-cbind(unlist(distoutput)[output1[[3]]],unlist(distoutput)[output1[[3]]+1])
  }
  if(Nstress==3){
    Stressset<-cbind(unlist(distoutput)[output1[[3]]],unlist(distoutput)[output1[[3]]+1],unlist(distoutput)[output1[[3]]+2])
  }
  setvect <- unlist(distoutput)[output1[[5]]]
  distparam0<-mean(setvect[!is.na(setvect)])

  if (ls=="Linear"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    lsoutput <- function(Lest){
      params  <- lm(Lest ~ poly(Stressset, 1, raw=TRUE))
      lsparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      (theta[2]+Sfrom*theta[1])/(theta[2]+Sto*theta[1])
    }
    adjparam <- function(theta) {
      col1<-replace(theta*StoAdj - SfromAdj,iblocked2,0)
      col2<-replace(theta - 1,iblocked2,0)
      mati<-matrix(c(col1[col1!=0],col2[col2!=0]),nrow=remNsteps,ncol=2,byrow=FALSE)
      v1<-pinv(mati)%*%rep(0,remNsteps)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="Exponential"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    lsoutput <- function(Lest){
      params  <- lm(log(Lest) ~ poly(Stressset, 1, raw=TRUE))
      lsparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      exp(theta[1]*(Sfrom-Sto))
    }
    adjparam <- function(theta) {
      v1<-log(theta)/(SfromAdj - StoAdj)
      v1<-replace(v1,iblocked2,0)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="Arrhenius"){
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b, lsparams[3] - R^2
    # Temperature HAS to be in Kelvin for this to work
    K<-8.617385e-5
    lsoutput <- function(Lest){
      params  <- lm(log(Lest) ~ poly(1/Stressset, 1, raw=TRUE))
      lsparams <- c(K*summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      exp((theta[1]/K)*((1/Sfrom)-(1/Sto)))
    }
    adjparam <- function(theta) {
      v1<-(K*log(theta))/((1/SfromAdj) - (1/StoAdj))
      v1<-replace(v1,iblocked2,0)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="Eyring"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    lsoutput <- function(Lest){
      Lvals<-log(Lest)
      params  <- nls(Lvals ~ log(b) -log(Stressset) + (a/Stressset),start = list(a = 1,b = 3))
      lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      (Sto/Sfrom)*exp(theta[1]*((1/Sfrom)-(1/Sto)))
    }
    adjparam <- function(theta) {
      v1<-log(theta*(SfromAdj/StoAdj))/((1/SfromAdj) - (1/StoAdj))
      v1<-replace(v1,iblocked2,0)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="Eyring2"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    lsoutput <- function(Lest){
      Lvals<-log(Lest)
      params  <- nls(Lvals ~ -log(Stressset) + (b/Stressset) - a,start = list(a = 1,b = 3))
      lsparams <- c(summary(params)$coefficients[1,1],summary(params)$coefficients[2,1])
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      (Sto/Sfrom)*exp(theta[2]*((1/Sfrom)-(1/Sto)))
    }
    adjparam <- function(theta) {
      v1<-log(theta*(SfromAdj/StoAdj))/((1/SfromAdj) - (1/StoAdj))
      v1<-replace(v1,iblocked2,0)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[2]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="Power"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    lsoutput <- function(Lest){
      params  <- lm(log(Lest) ~ poly(log(Stressset), 1, raw=TRUE))
      lsparams <- c(summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      (Sfrom/Sto)^lsparams[1]
    }
    adjparam <- function(theta) {
      v1<-log(theta)/log(SfromAdj/StoAdj)
      v1<-replace(v1,iblocked2,0)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="InversePower"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    #positivity_v[ishift+2]<-1

    lsoutput <- function(Lest){
      params  <- lm(log(Lest) ~ poly(log(Stressset), 1, raw=TRUE))
      lsparams <- c(-summary(params)$coefficients[2,1],exp(summary(params)$coefficients[1,1]))
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta,Sfrom,Sto) {
      (Sto/Sfrom)^theta[1]
    }
    adjparam <- function(theta,SfromAdj,StoAdj,iblocked2) {
      v1<-log(theta)/log(StoAdj/SfromAdj)
      v1<-replace(v1,iblocked2,0)
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-mean(theta)
      return(newlsparams)
    }
  }
  if (ls=="Logarithmic"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - R^2
    lsoutput <- function(Lest){
      params  <- lm(Lest ~ poly(log(Stressset), 1, raw=TRUE))
      lsparams <- c(summary(params)$coefficients[2,1],summary(params)$coefficients[1,1])
      R2 <- summary(params)$r.squared
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      (theta[2]+log(Sfrom)*theta[1])/(theta[2]+log(Sto)*theta[1])
    }
  }
  if (ls=="MultiStress"){
    # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
    if(Nstress<2) {
      stop('Select a data set with more than one stress type.')
    }
    lsoutput <- function(Lest){
      Lvals<-log(Lest)
      Svals<-matrix(c(rep(1,length(Stressset[,1])),Stressset),nrow=remNsteps,ncol=1+Nstress,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      lsparams <- c(params)
      R2 <- 1
      return(list(lsparams,R2))
    }
  }
  if (ls=="TempHumidity"){
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
    if(Nstress<2) {
      stop('Select a data set with more than one stress type.')
    }
    lsoutput <- function(Lest){
      Lvals<-log(Lest)
      Svals<-matrix(c(rep(1,length(Stressset[,1])),1/Stressset[,therm],1/Stressset[,alttherm]),nrow=remNsteps,ncol=3,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      params[1]<-exp(params[1])
      lsparams <- c(params)
      R2 <- 1
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      exp((theta[2]*((1/Sfrom[[1]]) - (1/Sto[[1]])))+(theta[3]*((1/Sfrom[[2]]) - (1/Sto[[2]]))))
    }
    adjparam <- function(theta) {
      b<-replace(theta,iblocked2,0)
      col1<-replace(SfromAdj[[1]] - StoAdj[[1]],iblocked2,0)
      col2<-replace(SfromAdj[[2]] - StoAdj[[2]],iblocked2,0)
      # mati<-matrix(c(col1[col1!=0],col2[col2!=0]),nrow=remNsteps,ncol=2,byrow=FALSE)
      mati<-matrix(c(col1[col1!=0],col2[col2!=0]),nrow=length(col1[col1!=0]),ncol=2,byrow=FALSE)
      return(list(b,col1,col2,mati))
      v1<-pinv(mati)%*%log(b[b!=0])
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[2]<-theta[1]
      newlsparams[3]<-theta[2]
      return(newlsparams)
    }
  }
  if (ls=="TempNonthermal"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    if(Nstress<2) {
      stop('Select a data set with more than one stress type.')
    }
    lsoutput <- function(Lest){
      Lvals<-log(Lest)
      Svals<-matrix(c(1/Stressset[,alttherm],-log(Stressset[,therm]),rep(1,length(Stressset[,1]))),nrow=remNsteps,ncol=3,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      params[3]<-exp(params[3])
      lsparams <- c(params)
      R2 <- 1
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      ((Sto[[1]]/Sfrom[[1]])^theta[2])*exp(-theta[1]*((1/Sto[[2]])-(1/Sfrom[[2]])))
    }
    adjparam <- function(theta) {
      b<-replace(theta,iblocked2,0)
      col1<-replace((1/SfromAdj[[2]]) - (1/StoAdj[[2]]),iblocked2,0)
      col2<-replace(log(StoAdj[[1]]) - log(SfromAdj[[1]]),iblocked2,0)
      mati<-matrix(c(col1[col1!=0],col2[col2!=0]),nrow=length(col1[col1!=0]),ncol=2,byrow=FALSE)
      v1<-pinv(mati)%*%log(b[b!=0])
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[1]<-theta[1]
      newlsparams[2]<-theta[2]
      return(newlsparams)
    }
  }
  if (ls=="Eyring3"){
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d
    if(Nstress<2) {
      stop('Select a data set with more than one stress type.')
    }
    lsoutput <- function(Lest){
      Lvals<-log(Lest)+log(Stressset[,therm])
      Svals<-matrix(c(rep(1,length(Stressset[,1])),1/Stressset[,therm],Stressset[,alttherm],Stressset[,alttherm]/Stressset[,therm]),nrow=remNsteps,ncol=4,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      lsparams <- c(params)
      R2 <- 1
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      (Sto[[1]]/Sfrom[[1]])*exp((theta[2]*((1/Sfrom[[1]])-(1/Sto[[1]])))+(theta[3]*(Sfrom[[2]]-Sto[[2]]))+(theta[4]*((Sfrom[[2]]/Sfrom[[1]])-(Sto[[2]]/Sto[[1]]))))
    }
    adjparam <- function(theta) {
      b<-replace(theta,iblocked2,0)
      c<-replace(StoAdj[[1]]/SfromAdj[[1]],iblocked2,0)
      col1<-replace((1/SfromAdj[[1]]) - (1/StoAdj[[1]]),iblocked2,0)
      col2<-replace(SfromAdj[[2]] - StoAdj[[2]],iblocked2,0)
      col3<-replace((SfromAdj[[2]]/SfromAdj[[1]]) - (StoAdj[[2]]/StoAdj[[1]]),iblocked2,0)
      mati<-matrix(c(col1[col1!=0],col2[col2!=0],col3[col3!=0]),nrow=length(col1[col1!=0]),ncol=3,byrow=FALSE)
      v1<-pinv(mati)%*%(log(b[b!=0])-log(c[c!=0]))
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[2]<-theta[1]
      newlsparams[3]<-theta[2]
      newlsparams[4]<-theta[3]
      return(newlsparams)
    }
  }

  if (ls=="PH1"){
    # lsparams[1] - parameter beta_0, lsparams[2] - parameter beta_1, lsparams[3] - beta_2, lsparams[4] - R^2
    lsoutput <- function(Lest){
      Lvals<-log(Lest)
      Svals<-matrix(c(rep(-1,length(Stressset[,1])),-1/Stressset[,therm],-1/Stressset[,alttherm]),nrow=remNsteps,ncol=3,byrow=FALSE)
      params  <- pinv(Svals)%*%Lvals
      lsparams <- c(params)
      R2 <- 1
      return(list(lsparams,R2))
    }
    init_AFn <- function(theta) {
      exp((theta[2]*((1/Sto[[1]]) - (1/Sfrom[[1]])))+(theta[3]*((1/Sto[[2]]) - (1/Sfrom[[2]]))))
      # exp((theta[2]*((1/Sfrom[[1]]) - (1/Sto[[1]])))+(theta[3]*((1/Sfrom[[2]]) - (1/Sto[[2]]))))
    }
    adjparam <- function(theta) {
      b<-replace(theta,iblocked2,0)
      col1<-replace(StoAdj[[1]] - SfromAdj[[1]],iblocked2,0)
      col2<-replace(StoAdj[[2]] - SfromAdj[[2]],iblocked2,0)
      # mati<-matrix(c(col1[col1!=0],col2[col2!=0]),nrow=remNsteps,ncol=2,byrow=FALSE)
      mati<-matrix(c(col1[col1!=0],col2[col2!=0]),nrow=length(col1[col1!=0]),ncol=2,byrow=FALSE)
      # return(list(b,col1,col2,mati))
      v1<-pinv(mati)%*%log(b[b!=0])
      v1<-v1[v1!=0]
      return(v1)
    }
    updateparam <- function(theta) {
      newlsparams<-lsparams
      newlsparams[2]<-theta[1]
      newlsparams[3]<-theta[2]
      return(newlsparams)
    }
  }

  # Setup the matrix blocks to update step-stress times
  if(is.list(data[[1]])==FALSE){
    Nsteps<-length(which(unique(cumdmg)>0))
    totNsteps<-dim(stepstresstable)[1]
    remNsteps<-dim(stpstrdatsort[[2]])[1]
    iblocked<-which(t(matrix(rep(1,totNsteps^2),totNsteps,totNsteps)*c(1:totNsteps))/(matrix(rep(1,totNsteps^2),totNsteps,totNsteps)*c(1:totNsteps))<=1)
    iblocked2<-which(t(matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*c(1:remNsteps))/(matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*c(1:remNsteps))<=1)
    stepsno<-as.numeric(rownames(stpstrdatsort[[2]]))

    # Initialize stress matrices Sfrom and Sto
    if(Nstress==1){
      Sfrom<-matrix(rep(1,totNsteps^2),totNsteps,totNsteps)*stepstresstable[,1]
      Sto<-t(Sfrom)
      tequiv<-matrix(rep(1,totNsteps^2),totNsteps,totNsteps)*stepstresstable[,dim(stepstresstable)[2]]
      # For Adjusting
      SfromAdj<-matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*stpstrdatsort[[2]][,1]
      StoAdj<-t(SfromAdj)
      cumdmgfrom<-matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*cumdmg[stepsno]
      cumdmgto<-t(cumdmgfrom)
    }
    if(Nstress>=2){
      Sfrom<-vector("list",Nstress)
      Sto<-vector("list",Nstress)
      tequiv<-matrix(rep(1,totNsteps^2),totNsteps,totNsteps)*stepstresstable[,dim(stepstresstable)[2]]
      # For Adjusting
      SfromAdj<-vector("list",Nstress)
      StoAdj<-vector("list",Nstress)
      cumdmgfrom<-matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*cumdmg[stepsno]
      cumdmgto<-t(cumdmgfrom)
      for(i2 in 1:Nstress){
        Sfrom[[i2]]<-matrix(rep(1,totNsteps^2),totNsteps,totNsteps)*stepstresstable[,i2]
        Sto[[i2]]<-t(Sfrom[[i2]])
        SfromAdj[[i2]]<-matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*stpstrdatsort[[2]][,i2]
        StoAdj[[i2]]<-t(SfromAdj[[i2]])
      }
    }
  } else{
    Nsteps<-vector("list",length(cumdmg))
    totNsteps<-vector("list",length(cumdmg))
    remNsteps<-vector("list",length(cumdmg))
    iblocked<-vector("list",length(cumdmg))
    iblocked2<-vector("list",length(cumdmg))
    stepsno<-vector("list",length(cumdmg))
    Sfrom<-vector("list",length(cumdmg))
    Sto<-vector("list",length(cumdmg))
    tequiv<-vector("list",length(cumdmg))
    SfromAdj<-vector("list",length(cumdmg))
    StoAdj<-vector("list",length(cumdmg))
    cumdmgfrom<-vector("list",length(cumdmg))
    cumdmgto<-vector("list",length(cumdmg))

    for(i2 in 1:length(cumdmg)){
      Nsteps[[i2]]<-length(which(unique(cumdmg[[i2]])>0))
      totNsteps[[i2]]<-dim(stepstresstable[[i2]])[1]
      remNsteps[[i2]]<-dim(stpstrdatsort[[i2]][[2]])[1]
      iblocked[[i2]]<-which(t(matrix(rep(1,totNsteps[[i2]]^2),totNsteps[[i2]],totNsteps[[i2]])*c(1:totNsteps[[i2]]))/(matrix(rep(1,totNsteps[[i2]]^2),totNsteps[[i2]],totNsteps[[i2]])*c(1:totNsteps[[i2]]))<=1)
      iblocked2[[i2]]<-which(t(matrix(rep(1,remNsteps[[i2]]^2),remNsteps[[i2]],remNsteps[[i2]])*c(1:remNsteps[[i2]]))/(matrix(rep(1,remNsteps[[i2]]^2),remNsteps[[i2]],remNsteps[[i2]])*c(1:remNsteps[[i2]]))<=1)
      stepsno[[i2]]<-as.numeric(rownames(stpstrdatsort[[i2]][[2]]))
      # Initialize stress matrices Sfrom and Sto
      if(Nstress==1){
        Sfrom[[i2]]<-matrix(rep(1,totNsteps[[i2]]^2),totNsteps[[i2]],totNsteps[[i2]])*stepstresstable[[i2]][,1]
        Sto[[i2]]<-t(Sfrom[[i2]])
        tequiv[[i2]]<-matrix(rep(1,totNsteps[[i2]]^2),totNsteps[[i2]],totNsteps[[i2]])*stepstresstable[[i2]][,dim(stepstresstable[[i2]])[2]]
        # For Adjusting
        SfromAdj[[i2]]<-matrix(rep(1,remNsteps[[i2]]^2),remNsteps[[i2]],remNsteps[[i2]])*stpstrdatsort[[i2]][[2]][,1]
        StoAdj[[i2]]<-t(SfromAdj[[i2]])
        cumdmgfrom[[i2]]<-matrix(rep(1,remNsteps[[i2]]^2),remNsteps[[i2]],remNsteps[[i2]])*cumdmg[[i2]][stepsno[[i2]]]
        cumdmgto[[i2]]<-t(cumdmgfrom[[i2]])
      }
    }
  }

  # Initialize life-stress model parameters
  lsparams0<-lsoutput(lifeest)[[1]]
  lsparams<-lsparams0
  lsparamsfirst<-lsparams0
  paramsfirst<-c(distparam0,lsparamsfirst)

  # Initialize first updated parameter set for optimization
  if(is.list(data[[1]])==FALSE){
    AFn<-replace(init_AFn(lsparams0),iblocked,0)
    Tequiv<-colSums(replace(tequiv/AFn,iblocked,0))+stepstresstable[,dim(stepstresstable)[2]]
    Tequivfrom<-matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*Tequiv[stepsno]
    Tequivto<-t(Tequivfrom)
    AFAdj<-(cumdmgto/cumdmgfrom)*(Tequivfrom/Tequivto)
    AFAdj<-replace(AFAdj,iblocked2,0)
    nAdj<-adjparam(AFAdj,SfromAdj,StoAdj,iblocked2)
  } else{
    AFn<-vector("list",length(cumdmg))
    Tequiv<-vector("list",length(cumdmg))
    Tequivfrom<-vector("list",length(cumdmg))
    Tequivto<-vector("list",length(cumdmg))
    AFAdj<-vector("list",length(cumdmg))
    AFAdj<-vector("list",length(cumdmg))
    nAdj<-vector("list",length(cumdmg))
    for(i2 in 1:length(cumdmg)){
      AFn[[i2]]<-replace(init_AFn(lsparams0,Sfrom[[i2]],Sto[[i2]]),iblocked[[i2]],0)
      Tequiv[[i2]]<-colSums(replace(tequiv[[i2]]/AFn[[i2]],iblocked[[i2]],0))+stepstresstable[[i2]][,dim(stepstresstable[[i2]])[2]]
      Tequivfrom[[i2]]<-matrix(rep(1,remNsteps[[i2]]^2),remNsteps[[i2]],remNsteps[[i2]])*Tequiv[[i2]][stepsno[[i2]]]
      Tequivto[[i2]]<-t(Tequivfrom[[i2]])
      AFAdj[[i2]]<-(cumdmgto[[i2]]/cumdmgfrom[[i2]])*(Tequivfrom[[i2]]/Tequivto[[i2]])
      AFAdj[[i2]]<-replace(AFAdj[[i2]],iblocked2[[i2]],0)
      nAdj[[i2]]<-adjparam(AFAdj[[i2]],SfromAdj[[i2]],StoAdj[[i2]],iblocked2[[i2]])
    }
  }

  for (i3 in 1:100){
    if(is.list(data[[1]])==FALSE){
      AFn<-replace(init_AFn(lsparams0),iblocked,0)
      Tequiv<-colSums(replace(tequiv/AFn,iblocked,0))+stepstresstable[,dim(stepstresstable)[2]]
      Tequivfrom<-matrix(rep(1,remNsteps^2),remNsteps,remNsteps)*Tequiv[stepsno]
      Tequivto<-t(Tequivfrom)
      AFAdj<-(cumdmgto/cumdmgfrom)*(Tequivfrom/Tequivto)
      AFAdj<-replace(AFAdj,iblocked2,0)
      nAdj<-adjparam(AFAdj,SfromAdj,StoAdj,iblocked2)
      lsparams0<-updateparam(nAdj)
    } else{
      AFn<-vector("list",length(cumdmg))
      Tequiv<-vector("list",length(cumdmg))
      Tequivfrom<-vector("list",length(cumdmg))
      Tequivto<-vector("list",length(cumdmg))
      AFAdj<-vector("list",length(cumdmg))
      AFAdj<-vector("list",length(cumdmg))
      nAdj<-vector("list",length(cumdmg))
      for(i2 in 1:length(cumdmg)){
        AFn[[i2]]<-replace(init_AFn(lsparams0,Sfrom[[i2]],Sto[[i2]]),iblocked[[i2]],0)
        Tequiv[[i2]]<-colSums(replace(tequiv[[i2]]/AFn[[i2]],iblocked[[i2]],0))+stepstresstable[[i2]][,dim(stepstresstable[[i2]])[2]]
        Tequivfrom[[i2]]<-matrix(rep(1,remNsteps[[i2]]^2),remNsteps[[i2]],remNsteps[[i2]])*Tequiv[[i2]][stepsno[[i2]]]
        Tequivto[[i2]]<-t(Tequivfrom[[i2]])
        AFAdj[[i2]]<-(cumdmgto[[i2]]/cumdmgfrom[[i2]])*(Tequivfrom[[i2]]/Tequivto[[i2]])
        AFAdj[[i2]]<-replace(AFAdj[[i2]],iblocked2[[i2]],0)
        nAdj[[i2]]<-adjparam(AFAdj[[i2]],SfromAdj[[i2]],StoAdj[[i2]],iblocked2[[i2]])
        lsparams0<-updateparam(unlist(nAdj))
      }
    }
  }

  if(is.list(data[[1]])==FALSE){
    del_T<-Tequiv-stepstresstable[,dim(stepstresstable)[2]]
    del_T_trim<-del_T[stepsno]
    del_Tvector<-rep(0,length(stpstrdatsort[[1]][,1]))
    updatedata<-stpstrdatsort[[1]]

    if(Nstress==1){
      for(i2 in 1:length(stpstrdatsort[[2]][,1])){
        del_Tvector[which(stpstrdatsort[[1]][,3]==stpstrdatsort[[2]][i2,1])]<-del_T_trim[i2]
      }
    }
    if(Nstress==2){
      for(i2 in 1:length(stpstrdatsort[[2]][,1])){
        del_Tvector[intersect(which(stpstrdatsort[[1]][,3]==stpstrdatsort[[2]][i2,1]),which(stpstrdatsort[[1]][,4]==stpstrdatsort[[2]][i2,2]))]<-del_T_trim[i2]
      }
    }
    if(Nstress>=3){
      for(i2 in 1:length(stpstrdatsort[[2]][,1])){
        stressmatch<-vector(mode = "list",Nstress)
        for(i3 in 1:Nstress){
          stressmatch[[i3]]<-which(stpstrdatsort[[1]][,i3+2]==stpstrdatsort[[2]][i2,i3])
        }
        del_Tvector[Reduce(intersect, stressmatch)]<-del_T_trim[i2]
      }
    }
    updatedata[,1]<-updatedata[,1]+del_Tvector
  } else{
    del_T<-vector("list",length(cumdmg))
    del_T_trim<-vector("list",length(cumdmg))
    del_Tvector<-vector("list",length(cumdmg))
    updatedata<-full_stpstrdata
    # updatestpstrdatsort<-vector("list",length(cumdmg))
    for(i3 in 1:length(cumdmg)){
      del_T[[i3]]<-Tequiv[[i3]]-stepstresstable[[i3]][,dim(stepstresstable[[i3]])[2]]
      del_T_trim[[i3]]<-del_T[[i3]][stepsno[[i3]]]
      del_Tvector[[i3]]<-rep(0,length(stpstrdatsort[[i3]][[1]][,1]))
      # updatedata[[i3]]<-stpstrdatsort[[i3]][[1]]
      if(Nstress==1){
        for(i2 in 1:length(stpstrdatsort[[i3]][[2]][,1])){
          del_Tvector[[i3]][which(stpstrdatsort[[i3]][[1]][,3]==stpstrdatsort[[i3]][[2]][i2,1])]<-del_T_trim[[i3]][i2]
        }
      }
      if(Nstress==2){
        for(i2 in 1:length(stpstrdatsort[[i3]][[2]][,1])){
          del_Tvector[[i3]][intersect(which(stpstrdatsort[[i3]][[1]][,3]==stpstrdatsort[[i3]][[2]][i2,1]),which(stpstrdatsort[[i3]][[1]][,4]==stpstrdatsort[[i3]][[2]][i2,2]))]<-del_T_trim[[i3]][i2]
        }
      }
      if(Nstress>=3){
        for(i2 in 1:length(stpstrdatsort[[i3]][[2]][,1])){
          stressmatch<-vector(mode = "list",Nstress)
          for(i4 in 1:Nstress){
            stressmatch[[i4]]<-which(stpstrdatsort[[i3]][[1]][,i4+2]==stpstrdatsort[[i3]][[2]][i2,i4])
          }
          del_Tvector[[i3]][Reduce(intersect, stressmatch)]<-del_T_trim[[i3]][i2]
        }
      }
      # updatedata[[i3]][,1]<-updatedata[[i3]][,1]+del_Tvector[[i3]]
      # updatestpstrdatsort[[i3]]<-stepstress.data(updatedata[[i3]],stepstresstable[[i3]])
      # if(i3 == 1){
      #   updatefull_stpstrdata <- updatestpstrdatsort[[i3]][[1]]
      # } else{
      #   updatefull_stpstrdata <- merge(updatefull_stpstrdata,updatestpstrdatsort[[i3]][[1]], all = TRUE, sort=FALSE)
      # }
    }
    updatedata[,1]<-updatedata[,1]+unlist(del_Tvector)
    # updatedata[[i3]][,1]<-updatedata[[i3]][,1]+del_Tvector[[i3]]
    # colnames(updatefull_stpstrdata) <- colnames(data[[1]], do.NULL = FALSE, prefix = "Obs.")
  }

  # Tabulate LSQ optimized estimates
  output2<-distlifeest(updatedata)
  distoutput<-output2[[1]]
  lifeest<-output2[[2]]
  setvect <- unlist(distoutput)[output2[[5]]]
  distparam<-mean(setvect[!is.na(setvect)])
  lsparams<-lsoutput(lifeest)[[1]]

  return(list(paramsfirst,c(distparam0,lsparams0),c(distparam,lsparams),lifeest,full_stpstrdata,updatedata))
  # return(list(c(distparam0,lsparams0),c(distparam,lsparams),lifeest,Sfrom,Sto,AFn,AFAdj,nAdj,tequiv,Tequiv,del_T,updatedata,c(distparam0,lsparamsfirst),updatedata))
}
