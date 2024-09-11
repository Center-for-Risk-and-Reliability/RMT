# Normal Probability Plot
# Developed by Dr. Reuel Smith, 2021-2024

probplot.nor <- function(data,pp,xlabel1="X",confid=0.95,stpstr_i=NULL,MLE_i=NULL) {
  library(ggplot2)

  # Legend colors
  col_legend <- c("red","blue","darkgreen","violet","aquamarine","orange","pink","darkblue","lightgreen","yellow","green")
  # Legend shapes
  shape_legend <- c(0:25)

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
  Pticks <- qnorm(Pticks1/100,mean=0,sd=1)
  Pticks1label <- c(0.1,0.2,0.3,"",0.5,"","","","",1,2,3,"",5,"","","","",10*c(1:9),95,99,99.9)
  fcB <- qnorm(c(.001,0.999),mean=0,sd=1)

  if (!is.null(dim(databystress))){
    # Single Stress data
    ixirc<-sort.xircdata(data)
    xiRFblock<-plotposit.select(ixirc[[2]],ixirc[[3]],pp)
    # Data pull
    XB<-xiRFblock[,1]
    FB <- qnorm(xiRFblock[,2],mean=0,sd=1)
    ttfc <- probplotparam.nor(xiRFblock[,1],xiRFblock[,2])
    ttfcrange <- ttfc[[1]]
    params<-ttfc[[3]]
    outputpp<-list(ttfc[[3]],ttfc[[4]])
  } else {
    ixirc_list<- vector(mode = "list", length = length(databystress))
    xiRFblock_list<- vector(mode = "list", length = length(databystress))
    XB_list<- vector(mode = "list", length = length(databystress))
    FB_list<- vector(mode = "list", length = length(databystress))
    Xbound_list<- vector(mode = "list", length = length(databystress))
    Xboundlow_list<- vector(mode = "list", length = length(databystress))
    Fboundlow_list<- vector(mode = "list", length = length(databystress))
    Xboundhigh_list<- vector(mode = "list", length = length(databystress))
    Fboundhigh_list<- vector(mode = "list", length = length(databystress))
    ttfc_list<- vector(mode = "list", length = length(databystress))
    ttfc_MLE_list<- vector(mode = "list", length = length(databystress))
    outputpp <- vector(mode = "list", length = 3*length(databystress)+1)
    SSEtot<-rep(1, length(databystress))
    ttfcrange<-c(1)
    XBfull<-c(1)
    FBfull<-c(1)
    # Multi-stress data
    for(i in 1:length(databystress)){
      ixirc_list[[i]]<-sort.xircdata(databystress[[i]])
      xiRFblock_list[[i]]<-plotposit.select(ixirc_list[[i]][[2]],ixirc_list[[i]][[3]],pp)
      # Data pull
      XB_list[[i]]<-xiRFblock_list[[i]][,1]
      FB_list[[i]] <- qnorm(xiRFblock_list[[i]][,2],mean=0,sd=1)
      ttfc_list[[i]]<-probplotparam.nor(xiRFblock_list[[i]][,1],xiRFblock_list[[i]][,2])
      if(is.null(MLE_i) == FALSE){
        # Compute MLE
        ttfc_MLE_list[[i]] <- distribution.MLEest(c(ttfc_list[[i]][[3]]),"Normal",ixirc_list[[i]][[2]],ixirc_list[[i]][[3]])
        # Recalculate x-points if MLE
        ttfc_list[[i]][[1]] <- c(ttfc_MLE_list[[i]][[1]][2])*fcB + c(ttfc_MLE_list[[i]][[1]][1])
      }
      outputpp[[i*3-2]]<-databystress[[i]][1,3:length(databystress[[i]][1,])]
      if(is.null(MLE_i) == TRUE){
        outputpp[[i*3-1]]<-ttfc_list[[i]][[3]]
      } else{
        outputpp[[i*3-1]]<-ttfc_MLE_list[[i]][[1]]
        ttfc_list[[i]][[3]]<-ttfc_MLE_list[[i]][[1]]
      }
      outputpp[[i*3]]<-ttfc_list[[i]][[4]]
      SSEtot[i]<-ttfc_list[[i]]$SSE
      XBfull<-c(XBfull,XB_list[[i]])
      FBfull<-c(FBfull,FB_list[[i]])
      ttfcrange<-c(ttfcrange,ttfc_list[[i]][[1]])
    }
    outputpp[[3*length(databystress)+1]]<-singledat
    ttfcrange<-ttfcrange[2:length(ttfcrange)]
    XBfull<-XBfull[2:length(XBfull)]
    FBfull<-FBfull[2:length(FBfull)]
  }

  # Calculate the upper and lower bounds
  Fbound <- c(linspace(0.0001,0.0009,9),linspace(.001,.999,982),linspace(0.9991,0.9999,9))
  FB_bound <- qnorm(Fbound,mean=0,sd=1)
  N <- dim(data)[1]
  Z <- qnorm(1-(1-confid)/2,0,1)
  Fdiff <- Z*sqrt((Fbound*(1-Fbound))/N)
  F_high <- Fbound + Fdiff
  F_high[which(F_high>1)] <- 0.9999999
  FB_high <- qnorm(F_high,mean=0,sd=1)
  F_low <- Fbound - Fdiff
  F_low[which(F_low<0)] <- 0.00000001
  FB_low <- qnorm(F_low,mean=0,sd=1)

  for(i in 1:length(databystress)){
    Xbound_list[[i]] <- c(ttfc_list[[i]][[3]][2])*FB_bound + c(ttfc_list[[i]][[3]][1])
    Fboundlow_list[[i]] <- FB_low
    Xboundlow_list[[i]] <- Xbound_list[[i]]
    Fboundhigh_list[[i]] <- FB_high
    Xboundhigh_list[[i]] <- Xbound_list[[i]]
  }

  # ========================================================
  # Plotting
  # ========================================================
  if (!is.null(dim(databystress))){
    # Single Stress
    df <- data.frame(XScale = XB, Fscale = FB, data = rep("data",length(XB)))
    df2 <- data.frame(Xline = ttfc[[1]], Fline = ttfc[[2]], best_fit = rep("Best-fit",2))

    plotout<-ggplot() +
      geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
      scale_x_continuous(limits = min(ttfcrange), max(ttfcrange), breaks=Pticks1X, labels=Pticks1Xlabel) +
      scale_y_continuous(limits = c(min(fcB), max(fcB)), breaks=Pticks, labels=Pticks1label) +
      xlab(xlabel1) +
      ylab("Percent Failure")
    plotout <- plotout + geom_line(data=df2, aes(Xline,Fline, colour = best_fit), size = 0.9, linetype = "dashed")

  } else {
    # Multi-Stress
    data_legend <- logical(0)
    xlines <- rep(0,length(databystress)*2)
    Flines <- rep(0,length(databystress)*2)
    xlinesconf <- rep(0,length(databystress)*2001)
    xlinesconf_up <- rep(0,length(databystress)*1000)
    xlinesconf_down <- rep(0,length(databystress)*1000)
    Flinesconf <- rep(0,length(databystress)*2001)
    Flinesconf_up <- rep(0,length(databystress)*1000)
    Flinesconf_down <- rep(0,length(databystress)*1000)
    line_legend <- rep(0,length(databystress)*2)
    conf_legend <- rep(0,length(databystress)*2001)

    for(i in 1:length(databystress)){
      xlines[((i*2) - 1):(i*2)] <- ttfc_list[[i]][[1]]
      Flines[((i*2) - 1):(i*2)] <- ttfc_list[[i]][[2]]
      xlinesconf[((i*2001) - 2000):(i*2001)] <- c(Xboundlow_list[[i]],NA,Xboundhigh_list[[i]])
      Flinesconf[((i*2001) - 2000):(i*2001)] <- c(Fboundlow_list[[i]],NA,Fboundhigh_list[[i]])
      Flinesconf_up[((i*1000) - 999):(i*1000)] <- Fboundhigh_list[[i]]
      Flinesconf_down[((i*1000) - 999):(i*1000)] <- Fboundlow_list[[i]]
      xlinesconf_down[((i*1000) - 999):(i*1000)] <- Xboundlow_list[[i]]
      xlinesconf_up[((i*1000) - 999):(i*1000)] <- Xboundhigh_list[[i]]
      data_legend<-c(data_legend,rep(paste(c("Data for stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),length(XB_list[[i]])))
      line_legend[((i*2) - 1):(i*2)] <- rep(paste(c("Best-fit for stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),2)
      conf_legend[((i*2001) - 2000):(i*2001)] <- rep(paste(c(confid*100,"% Confidence for stress level ",databystress[[i]][1,3:length(databystress[[i]][1,])]),collapse = " "),2001)
    }
    df <- data.frame(XScale = XBfull, Fscale = FBfull, data = data_legend)
    df2 <- data.frame(Xline = xlines, Fline = Flines, best_fit = line_legend)
    df3 <- data.frame(Xline2up = xlinesconf_up, Xline2down = xlinesconf_down, Fline2up = Flinesconf_up, Fline2down = Flinesconf_down)

    # If step-stress problem, identify the index of comparison (RCS 01082024)
    if(is.null(stpstr_i)==FALSE){
      df <- df[which(stpstr_i==1),]
    }

    plotout<-ggplot() +
      geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
      scale_shape_manual(values=shape_legend[1:length(databystress)]) +
      scale_color_manual(values=col_legend[1:length(databystress)])+
      scale_x_continuous(limits = c(min(ttfcrange), max(ttfcrange))) +
      scale_y_continuous(limits = c(min(fcB), max(fcB)), breaks=Pticks, labels=Pticks1label) +
      xlab(xlabel1) +
      ylab("Percent Failure")

    plotout <- plotout + geom_line(data=df2, aes(Xline,Fline, colour = best_fit), size = 0.9, linetype = "dashed") +
      geom_ribbon(data=df3, aes(x=Xline2up, ymin=Fline2down, ymax=Fline2up), alpha=0.25, fill = "blue")
  }

  return(list(outputpp, SSEbyStress = SSEtot, summary.nonparametric=xiRFblock_list, prob_plot = plotout))
}
