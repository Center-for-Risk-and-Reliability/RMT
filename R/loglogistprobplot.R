# Log-logistic Probability Plot
# Developed by Dr. Reuel Smith, 2023

probplot.loglogist <- function(data,pp,xlabel1="X") {
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
  Pticks <- -log((1/(Pticks1/100))-1)
  Pticks1label <- c(0.1,0.2,0.3,"",0.5,"","","","",1,2,3,"",5,"","","","",10*c(1:9),95,99,99.9)
  fcB <- -log((1/c(.001,0.999))-1)

  if (!is.null(dim(databystress))){
    # Single Stress data
    ixirc<-sort.xircdata(data)
    xiRFblock<-plotposit.select(ixirc[[2]],ixirc[[3]],pp)
    # Data pull
    XB<-xiRFblock[,1]
    FB <- -log((1/xiRFblock[,2])-1)
    ttfc <- probplotparam.loglogist(xiRFblock[,1],xiRFblock[,2])
    ttfcrange <- ttfc[[1]]
    params<-ttfc[[3]]
    outputpp<-list(ttfc[[3]],ttfc[[4]])
  } else {
    ixirc_list<- vector(mode = "list", length = length(databystress))
    xiRFblock_list<- vector(mode = "list", length = length(databystress))
    XB_list<- vector(mode = "list", length = length(databystress))
    FB_list<- vector(mode = "list", length = length(databystress))
    ttfc_list<- vector(mode = "list", length = length(databystress))
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
      FB_list[[i]] <- -log((1/xiRFblock_list[[i]][,2])-1)
      ttfc_list[[i]]<-probplotparam.loglogist(xiRFblock_list[[i]][,1],xiRFblock_list[[i]][,2])
      outputpp[[i*3-2]]<-databystress[[i]][1,3:length(databystress[[i]][1,])]
      outputpp[[i*3-1]]<-ttfc_list[[i]][[3]]
      outputpp[[i*3]]<-ttfc_list[[i]][[4]]
      XBfull<-c(XBfull,XB_list[[i]])
      FBfull<-c(FBfull,FB_list[[i]])
      ttfcrange<-c(ttfcrange,ttfc_list[[i]][[1]])
    }
    outputpp[[3*length(databystress)+1]]<-singledat
    ttfcrange<-ttfcrange[2:length(ttfcrange)]
    XBfull<-XBfull[2:length(XBfull)]
    FBfull<-FBfull[2:length(FBfull)]
  }


  # Computes the upper and lower bound for the TTF axis in terms of log-time
  signs1 <- c(floor(log10(min(ttfcrange))):ceiling(log10(max(ttfcrange))))
  logtimes1 <- 10^signs1
  Pticks1X <- c(1:(9*length(logtimes1)-8))
  Pticks1X[1] <- logtimes1[1]
  Pticks1Xlabel <- Pticks1X

  for(i2 in 1:(length(signs1)-1)){
    Pticks1X[(9*i2-7):(9*(i2+1)-8)] <- logtimes1[i2]*c(2:10)
    Pticks1Xlabel[(9*i2-7):(9*(i2+1)-8)] <- c("","","","","","","","",logtimes1[i2+1])
  }


  # ========================================================
  # New Plotting
  # ========================================================
  if (!is.null(dim(databystress))){
    # Single Stress
    df <- data.frame(XScale = XB, Fscale = FB, data = rep("data",length(XB)))
    df2 <- data.frame(Xline = ttfc[[1]], Fline = -ttfc[[2]], best_fit = rep("Best-fit",2))

    plotout<-ggplot() +
      geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
      scale_x_continuous(limits = c(10^min(signs1), 10^max(signs1)), breaks=Pticks1X, labels=Pticks1Xlabel) +
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

    for(i in 1:length(databystress)){
      xlines[((i*2) - 1):(i*2)] <- ttfc_list[[i]][[1]]
      Flines[((i*2) - 1):(i*2)] <- -ttfc_list[[i]][[2]]
      data_legend<-c(data_legend,rep(paste(c("Data for stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),length(XB_list[[i]])))
      line_legend[((i*2) - 1):(i*2)] <- rep(paste(c("Best-fit for stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),2)
    }
    df <- data.frame(XScale = XBfull, Fscale = FBfull, data = data_legend)
    df2 <- data.frame(Xline = xlines, Fline = Flines, best_fit = line_legend)

    plotout<-ggplot() +
      geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
      scale_x_continuous(trans = 'log10', limits = c(10^min(signs1), 10^max(signs1)), breaks=Pticks1X, labels=Pticks1Xlabel) +
      scale_y_continuous(limits = c(min(fcB), max(fcB)), breaks=Pticks, labels=Pticks1label) +
      xlab(xlabel1) +
      ylab("Percent Failure")

    plotout <- plotout + geom_line(data=df2, aes(Xline,Fline, colour = best_fit), size = 0.9, linetype = "dashed")
  }

  # return(plotout)
  return(list(outputpp, prob_plot = plotout))
}
