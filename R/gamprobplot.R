# Gamma Probability Plot
# Developed by Dr. Reuel Smith, 2022-2023

probplot.gam <- function(data,pp,xlabel1="X") {
  library(pracma)
  library(zipfR)
  library(ggplot2)

  dcount<-checkdatacount(checkstress(data))
  # Error message to check if there is any plots for a life distribution estimate
  if(sum(dcount[[1]])==0){
    stop('Please check that there are at least two failure data that occur in the same stress level.')
  }

  # Check for whether the curve vector and censor vector are the same or not.  If
  # they are identical it means all relevant curves have failure data.  If not
  # then it is an indicator of some data not having curves and therefore must be
  # held over for further analysis.
  if(identical(dcount[[1]],dcount[[2]])){
    databystress<-checkstress(data)
    singledat<-NULL
  } else{
    databystress<-checkstress(data)[which(dcount[[1]] == 1 & dcount[[2]] == 1)]
    singledat<-checkstress(data)[which(!dcount[[1]] == 1 & dcount[[2]] == 1)]
  }

  # Initialize Percent Failure ticks
  Pticks1 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,c(1:10),10*c(2:9),95,99,99.9)
  Pticks1label <- c(0.1,0.2,0.3,"",0.5,"","","","",1,2,3,"",5,"","","","",10*c(1:9),95,99,99.9)

  # if (!is.null(dim(databystress))){
  if (is.list(databystress) && length(databystress)==1){
    # Single Stress data
    ixirc<-sort.xircdata(data)
    xiRFblock<-plotposit.select(ixirc[[2]],ixirc[[3]],pp)
    # Data pull
    XB<-xiRFblock[,1]

    ttfc <- probplotparam.gam(xiRFblock[,1],xiRFblock[,2])
    ttfcrange <- ttfc[[1]]
    params<-ttfc[[3]]
    FB <- Rgamma.inv(params[1],xiRFblock[,2])
    outputpp<-list(ttfc[[3]],ttfc[[4]])
    # These are now to be dependent on the window
    Pticks <- c(Rgamma.inv(params[1],Pticks1[which(Pticks1 < 90)]/100),Rgamma.inv(params[1],1-Pticks1[which(Pticks1 >= 90)]/100,lower = FALSE))
    fcB <- c(Rgamma.inv(params[1],0.001),Rgamma.inv(params[1],1-0.999,lower = FALSE))
  } else {
    ixirc_list<- vector(mode = "list", length = length(databystress))
    xiRFblock_list<- vector(mode = "list", length = length(databystress))
    XB_list<- vector(mode = "list", length = length(databystress))
    FB_list<- vector(mode = "list", length = length(databystress))
    ttfc_list<- vector(mode = "list", length = length(databystress))
    Pticks_list<- vector(mode = "list", length = length(databystress))
    fcB_list<- vector(mode = "list", length = length(databystress))
    outputpp <- vector(mode = "list", length = 3*length(databystress)+1)
    ttfcrange<-c(1)
    XBfull<-c(1)
    FBfull<-c(1)
    # Multi-stress data
    for(i in 1:length(databystress)){
      ixirc_list[[i]]<-sort.xircdata(databystress[[i]])
      xiRFblock_list[[i]]<-plotposit.select(ixirc_list[[i]][[2]],ixirc_list[[i]][[3]],pp)
      # Data pull
      XB_list[[i]]<-xiRFblock_list[[i]][,1]

      ttfc_list[[i]]<-probplotparam.gam(xiRFblock_list[[i]][,1],xiRFblock_list[[i]][,2])
      outputpp[[i*3-2]]<-databystress[[i]][1,3:length(databystress[[i]][1,])]
      outputpp[[i*3-1]]<-ttfc_list[[i]][[3]]
      outputpp[[i*3]]<-ttfc_list[[i]][[4]]
      FB_list[[i]] <- Rgamma.inv(ttfc_list[[i]][[3]][1],xiRFblock_list[[i]][,2])
      XBfull<-c(XBfull,XB_list[[i]])
      FBfull<-c(FBfull,FB_list[[i]])
      ttfcrange<-c(ttfcrange,ttfc_list[[i]][[1]])
      # These are now to be dependent on the window
      Pticks_list[[i]] <- Rgamma.inv(ttfc_list[[i]][[3]][1],Pticks1/100)
      fcB_list[[i]] <- c(Rgamma.inv(ttfc_list[[i]][[3]][1],0.001),Rgamma.inv(ttfc_list[[i]][[3]][1],1-.999,lower = FALSE))
    }
    outputpp[[3*length(databystress)+1]]<-singledat
    ttfcrange<-ttfcrange[2:length(ttfcrange)]
    XBfull<-XBfull[2:length(XBfull)]
    FBfull<-FBfull[2:length(FBfull)]
  }

  # ========================================================
  # Plotting
  # ========================================================
  # if (!is.null(dim(databystress))){
  if (is.list(databystress) && length(databystress)==1){
    # Single Stress
    df <- data.frame(XScale = XB, Fscale = FB, data = rep("data",length(XB)))
    df2 <- data.frame(Xline = ttfc[[1]], Fline = ttfc[[2]], best_fit = rep("Best-fit",2))

    plotout<-ggplot() +
      geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
      scale_y_continuous(limits = c(min(fcB), max(fcB)), breaks=Pticks, labels=Pticks1label) +
      xlab(xlabel1) +
      ylab("Percent Failure")
    plotout <- plotout + geom_line(data=df2, aes(Xline,Fline, colour = best_fit), size = 0.9, linetype = "dashed")

  } else {
    # Multi-Stress
    data_legend <- logical(0)
    xlines <- rep(0,length(databystress)*2)
    Flines <- rep(0,length(databystress)*2)
    line_legend <- rep(0,length(databystress)*2)
    plotout <- vector(mode = "list", length = length(databystress))

    for(i in 1:length(databystress)){
      # Generate a list of plots per stress level
      df <- data.frame(XScale = XB_list[[i]], Fscale = FB_list[[i]], data = rep(paste(c("Stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),length(XB_list[[i]])))
      df2 <- data.frame(Xline = ttfc_list[[i]][[1]], Fline = ttfc_list[[i]][[2]], best_fit = rep(paste(c("Stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),2))

      plotout[[i]]<-ggplot() +
        geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
        scale_y_continuous(limits = c(min(fcB_list[[i]]), max(fcB_list[[i]])), breaks=Pticks_list[[i]], labels=Pticks1label) +
        xlab(xlabel1) +
        ylab("Percent Failure")
      plotout[[i]] <- plotout[[i]] + geom_line(data=df2, aes(Xline,Fline, colour = best_fit), size = 0.9, linetype = "dashed")
    }
  }

  return(list(outputpp, prob_plot = plotout))
}
