# Exponential Probability Plot
# Developed by Dr. Reuel Smith, 2021-2025

probplot.exp <- function(data,pp="Blom",xlabel1="X",confid=0.95,
                         nobounds = NULL,stpstr_i=NULL,MLE_i=NULL,CDFdata = NULL,CDFrangesetting = 1,
                         stressunit1 = NULL, stressunit2 = NULL) {
  library(pracma)
  library(ggplot2)

  # Legend colors
  col_legend <- c("red","blue","darkgreen","violet","gold","orange","pink2","darkblue","lightgreen","yellow","green","darkviolet","darkorange","darkred","purple","royalblue","brown","lightpink","tan","darkgray","aquamarine","sienna","limegreen","mediumpurple3","chocolate","red4")
  # Legend shapes
  shape_legend <- c(0:25)
  # Legend line type
  linetype_legend <- rep(c(1,2,4,5,6),5)

  # State label for best fit line legend
  if(is.null(MLE_i)==TRUE){ # Switch off for LSQ
    best_fit_label <- "Exponential (LSQ) Best-Fit"
  }
  if(is.null(MLE_i)==FALSE){ # Switch on for MLE
    best_fit_label <- "Exponential (MLE) Best-Fit"
  }

  # 9/25/2025 - Check if data is a vector or not,  If so, augment the data such that it can be processed as
  # if it has the three column minimum needed to compute
  if(min(size(data)) == 1){
    data <- cbind(sort(data),rep(1,length(data)),rep(1,length(data)))
  }
  if(size(data)[2] == 2){ # add column if only the data and censored column are entered
    data <- cbind(data,rep(1,length(data[,1])))
  }

  dcount<-checkdatacount(checkstress(data))
  # Error message to check if there is any plots for a life distribution estimate
  # if(sum(dcount[[1]])==0){
  #   stop('Please check that there are at least two failure data that occur in the same stress level.')
  # }

  # Check for whether the curve vector and censor vector are the same or not.  If
  # they are identical it means all relevant curves have failure data.  If not
  # then it is an indicator of some data not having curves and therefore must be
  # held over for further analysis.
  if(identical(dcount[[1]],dcount[[2]])){
    databystress<-checkstress(data)
    singledat<-NULL
    if(is.null(CDFdata)==FALSE){
      CDFbystress <- CDFdata
    }
  } else{
    databystress<-checkstress(data)[which(dcount[[1]] == 1 & dcount[[2]] == 1)]
    singledat<-checkstress(data)[which(!dcount[[1]] == 1 & dcount[[2]] == 1)]
    if(is.null(CDFdata)==FALSE){
      CDFbystress <- CDFdata[which(!dcount[[1]] == 1 & dcount[[2]] == 1)]
    }
  }

  # Initialize Percent Failure ticks
  if(CDFrangesetting == 1){ # Minitab range 1% to 99%
    Fmin <- 0.01
    Fmax <- 0.99
    Pticks1 <- c(c(1:10),10*c(2:9),95,99)
    Pticks <- -log(1-Pticks1*0.01)
    Pticks1label <- c(1,"","","","","","","","",10*c(1:9),95,99)
    fcB <- -log(1-c(0.01,0.99))
  }
  if(CDFrangesetting == 2){ # Weibull++ range 0.1% to 99.9%
    Fmin <- 0.001
    Fmax <- 0.999
    Pticks1 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,c(1:10),10*c(2:9),95,99,99.9)
    Pticks <- -log(1-Pticks1*0.01)
    Pticks1label <- c(0.1,"","","","","","","","",1,"","","","","","","","",10*c(1:9),95,99,99.9)
    fcB <- -log(1-c(.001,0.999))
  }

  if (!is.null(dim(databystress))){
    # Single Stress data
    ixirc<-sort.xircdata(data)
    xiRFblock<-plotposit.select(ixirc[[2]],ixirc[[3]],pp)
    # Data pull
    XB<-xiRFblock[,1]
    FB <- -log(xiRFblock[,3])
    ttfc <- probplotparam.exp(xiRFblock[,1],xiRFblock[,2])
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
    outputpp1 <- vector(mode = "list", length = length(databystress))
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
      FB_list[[i]] <- -log(xiRFblock_list[[i]][,3])
      ttfc_list[[i]]<-probplotparam.exp(xiRFblock_list[[i]][,1],xiRFblock_list[[i]][,2],CDFrangesetting)
      if(is.null(MLE_i) == TRUE){  # Default for LSQ
        outputpp[[i*3-1]]<-ttfc_list[[i]][[3]] # LSQ parameters
        outputpp1[[i]][[2]]<-ttfc_list[[i]][[3]] # LSQ parameters
      }
      if(is.null(MLE_i) == FALSE){  # MLE for λ
        # Compute MLE
        ttfc_MLE_list[[i]] <- distribution.MLEest(c(ttfc_list[[i]][[3]]),"Exponential",ixirc_list[[i]][[2]],ixirc_list[[i]][[3]])
        # Recalculate x-points if MLE
        ttfc_list[[i]][[1]] <- fcB/ttfc_MLE_list[[i]][[1]][1]

        outputpp[[i*3-1]]<-matrix(ttfc_MLE_list[[i]][[1]], nrow = 1, ncol = 1, byrow = TRUE,dimnames = list(c("Exponential Parameter"),c("lambda")))
        outputpp1[[i]][[2]]<-matrix(ttfc_MLE_list[[i]][[1]], nrow = 1, ncol = 1, byrow = TRUE,dimnames = list(c("Exponential Parameter"),c("lambda")))
        ttfc_list[[i]][[3]]<-ttfc_MLE_list[[i]][[1]]
      }
      outputpp[[i*3-2]]<-databystress[[i]][1,3:length(databystress[[i]][1,])] # Stress level
      outputpp[[i*3]]<-ttfc_list[[i]][[4]] # Coefficient of determination
      SSEtot[i]<-ttfc_list[[i]]$SSE
      XBfull<-c(XBfull,XB_list[[i]])
      FBfull<-c(FBfull,FB_list[[i]])
      ttfcrange<-c(ttfcrange,ttfc_list[[i]][[1]])

      # Printed Output
      outputpp1[[i]][[1]]<-databystress[[i]][1,3:length(databystress[[i]][1,])] # Stress level
      if(is.null(MLE_i) == TRUE){
        outputpp1[[i]][[3]]<-ttfc_list[[i]]$SSE
        outputpp1[[i]][[4]]<-ttfc_list[[i]][[4]] # Coefficient of determination
        names(outputpp1[[i]]) <- c("Stress Level","Parameter Estimates","SSE","Coefficient of Determination (R2)")
      } else{
        outputpp1[[i]][[3]]<-matrix(unlist(ttfc_MLE_list[[i]][[3]]), nrow = 2, ncol = 1, byrow = FALSE,dimnames = list(c("Lower CI","Upper CI"),c("lambda"))) # Confidence
        outputpp1[[i]][[4]]<-ttfc_list[[i]]$SSE # SSE
        outputpp1[[i]][[5]]<-ttfc_MLE_list[[i]][[4]] # loglikelihood
        outputpp1[[i]][[6]]<-ttfc_MLE_list[[i]][[5]] # likelihood
        outputpp1[[i]][[7]]<-ttfc_MLE_list[[i]][[6]] # AIC
        outputpp1[[i]][[8]]<-ttfc_MLE_list[[i]][[7]] # BIC
        names(outputpp1[[i]]) <- c("Stress Level","Parameter Estimates","CI","SSE","loglikelihood","likelihood","AIC","BIC")
      }
    }

    outputpp[[3*length(databystress)+1]]<-singledat
    ttfcrange<-ttfcrange[2:length(ttfcrange)]
    XBfull<-XBfull[2:length(XBfull)]
    FBfull<-FBfull[2:length(FBfull)]
  }
  # Calculate the upper and lower bounds
  Z <- qnorm(1-(1-confid)/2,0,1)

  for(i in 1:length(databystress)){
    N <- dim(databystress[[i]])[1] # Compute N for stress level i
    SSElow <- function(F){
      (Fmin - F - Z*sqrt((F*(1 - F))/N))^2
    }
    SSEhigh <- function(F){
      (Fmax - F + Z*sqrt((F*(1 - F))/N))^2
    }
    # different for each stress level at least outside of the 982 points in the middle.  The upper and lower parts will
    # have their own value based on approximation
    A <- nlminb(Fmin,SSElow,hessian=TRUE,lower = 0,upper = 1)$par
    B <- nlminb(Fmax,SSEhigh,hessian=TRUE,lower = 0,upper = 1)$par

    # Draft Fbound for stress level i
    Fbound <- c(linspace(A,Fmin,250),linspace(Fmin,Fmax,500),linspace(Fmax,B,250))
    Fdiff <- Z*sqrt((Fbound*(1-Fbound))/N)
    F_high <- Fbound + Fdiff
    F_high[which(F_high<Fmin)] <- Fmin
    F_high[which(F_high>Fmax)] <- Fmax
    F_low <- Fbound - Fdiff
    F_low[which(F_low<Fmin)] <- Fmin
    F_low[which(F_low>Fmax)] <- Fmax

    FB_bound <- -log(1-Fbound)
    FB_high <- -log(1-F_high)
    FB_low <- -log(1-F_low)

    # F(x) = 1 - exp(-lam*x) => -lam*x = ln(1-F) => x = -ln(1-F)/lam
    Xbound_list[[i]] <- FB_bound/ttfc_list[[i]][[3]][1]
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
    data_breaks <- rep(0,length(databystress))
    xlines <- rep(0,length(databystress)*2)
    Flines <- rep(0,length(databystress)*2)
    xlinesconf <- rep(0,length(databystress)*2001)
    xlinesconf_up <- rep(0,length(databystress)*1000)
    xlinesconf_down <- rep(0,length(databystress)*1000)
    Flinesconf <- rep(0,length(databystress)*2001)
    Flinesconf_up <- rep(0,length(databystress)*1000)
    Flinesconf_down <- rep(0,length(databystress)*1000)
    line_legend <- rep(0,length(databystress)*2)
    line_breaks <- rep(0,length(databystress))
    conf_legend <- rep(0,length(databystress)*1000)
    conf_breaks <- rep(0,length(databystress))

    for(i in 1:length(databystress)){
      xlines[((i*2) - 1):(i*2)] <- ttfc_list[[i]][[1]]
      Flines[((i*2) - 1):(i*2)] <- ttfc_list[[i]][[2]]
      xlinesconf[((i*2001) - 2000):(i*2001)] <- c(Xboundlow_list[[i]],NA,Xboundhigh_list[[i]])
      Flinesconf[((i*2001) - 2000):(i*2001)] <- c(Fboundlow_list[[i]],NA,Fboundhigh_list[[i]])
      Flinesconf_up[((i*1000) - 999):(i*1000)] <- Fboundhigh_list[[i]]
      Flinesconf_down[((i*1000) - 999):(i*1000)] <- Fboundlow_list[[i]]
      xlinesconf_down[((i*1000) - 999):(i*1000)] <- Xboundlow_list[[i]]
      xlinesconf_up[((i*1000) - 999):(i*1000)] <- Xboundhigh_list[[i]]
      if(length(3:length(databystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        if(length(databystress) > 1){ # multiple groups of data
          data_legend<-c(data_legend,rep(paste(c("Data for",databystress[[i]][1,3:length(databystress[[i]][1,])]," units"),collapse = " "),length(XB_list[[i]])))
          data_breaks[i] <- paste(c("Data for ",databystress[[i]][1,3:length(databystress[[i]][1,])]," units"),collapse = " ")
          conf_legend[((i*1000) - 999):(i*1000)] <- rep(paste(c("% CI for ",databystress[[i]][1,3:length(databystress[[i]][1,])]," units"),collapse = " "),1000)
          conf_breaks[i] <- paste(c("% CI for ",databystress[[i]][1,3:length(databystress[[i]][1,])]," units"),collapse = " ")
        } else{ # only one group
          data_legend<-c(data_legend,rep("Data",length(XB_list[[i]])))
          data_breaks <- waiver()
          conf_legend[((i*1000) - 999):(i*1000)] <- rep("CI",1000)
          conf_breaks <- waiver()
        }
      }
      if(length(3:length(databystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,rep(paste(c("Data for",databystress[[i]][1,3:length(databystress[[i]][1,])]," ",stressunit1),collapse = " "),length(XB_list[[i]])))
        data_breaks[i] <- paste(c("Data for",databystress[[i]][1,3:length(databystress[[i]][1,])]," ",stressunit1),collapse = " ")
        conf_legend[((i*1000) - 999):(i*1000)] <- rep(paste(c("CI for",databystress[[i]][1,3:length(databystress[[i]][1,])],"",stressunit1),collapse = " "),1000)
        conf_breaks[i] <- paste(c("CI for",databystress[[i]][1,3:length(databystress[[i]][1,])],"",stressunit1),collapse = " ")
      }
      if(length(3:length(databystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for ",databystress[[i]][1,3]," units/",databystress[[i]][1,4]," units"),collapse = " "),length(XB_list[[i]])))
        data_breaks[i] <- paste(c("Data for ",databystress[[i]][1,3]," units/",databystress[[i]][1,4]," units"),collapse = " ")
        conf_legend[((i*1000) - 999):(i*1000)] <- rep(paste(c("CI for ",databystress[[i]][1,3]," units/",databystress[[i]][1,4]," units"),collapse = " "),1000)
        conf_breaks[i] <- paste(c("CI for ",databystress[[i]][1,3]," units/",databystress[[i]][1,4]," units"),collapse = " ")
      }
      if(length(3:length(databystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4],"",stressunit2),collapse = " "),length(XB_list[[i]])))
        data_breaks[i] <- paste(c("Data for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4],"",stressunit2),collapse = " ")
        conf_legend[((i*1000) - 999):(i*1000)] <- rep(paste(c("CI for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4],"",stressunit2),collapse = " "),1000)
        conf_breaks[i] <- paste(c("CI for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4],"",stressunit2),collapse = " ")
      }
      if(length(3:length(databystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4]," units"),collapse = " "),length(XB_list[[i]])))
        data_breaks[i] <- paste(c("Data for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4]," units"),collapse = " ")
        conf_legend[((i*1000) - 999):(i*1000)] <- rep(paste(c("CI for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4]," units"),collapse = " "),1000)
        conf_breaks[i] <- paste(c("CI for",databystress[[i]][1,3],"",stressunit1,"/",databystress[[i]][1,4]," units"),collapse = " ")
      }
      if(length(3:length(databystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for",databystress[[i]][1,3],""," units","/",databystress[[i]][1,4],"",stressunit2),collapse = " "),length(XB_list[[i]])))
        data_breaks[i] <- paste(c("Data for",databystress[[i]][1,3]," units","/",databystress[[i]][1,4],"",stressunit2),collapse = " ")
        conf_legend[((i*1000) - 999):(i*1000)] <- rep(paste(c("CI for",databystress[[i]][1,3]," units","/",databystress[[i]][1,4],"",stressunit2),collapse = " "),1000)
        conf_breaks[i] <- paste(c("CI for",databystress[[i]][1,3]," units","/",databystress[[i]][1,4],"",stressunit2),collapse = " ")
      }
      line_legend[((i*2) - 1):(i*2)] <- rep(paste(c("\U03BB = ",outputpp1[[i]][[2]],2),collapse = ""),2)
      if(length(databystress) > 1){ # multiple groups of data
        line_breaks[i] <- paste(c("\U03BB = ",outputpp1[[i]][[2]],2),collapse = "")
      } else{ # only one group
        line_breaks <- waiver()
      }
    }
    df <- data.frame(XScale = XBfull, Fscale = FBfull, data = data_legend)
    df2 <- data.frame(Xline = xlines, Fline = Flines, best_fit = line_legend)
    df3 <- data.frame(Xline2up = xlinesconf_up, Xline2down = xlinesconf_down, Fline2up = Flinesconf_up, Fline2down = Flinesconf_down,confidence = conf_legend)
    # return(list(df2,line_breaks))
    # If step-stress problem, identify the index of comparison (RCS 01082024)
    if(is.null(stpstr_i)==FALSE){
      df <- df[which(stpstr_i==1),]
    }

    plotout<-ggplot() +
      geom_point(data=df, aes(XScale,Fscale, shape = data), colour = 'black', size = 2.2) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_shape_manual("Raw Data",values=shape_legend[1:length(databystress)], breaks = data_breaks) +
      scale_x_continuous(expand=c(0, 0),limits = c(min(ttfcrange), max(ttfcrange))) +
      scale_y_continuous(expand=c(0, 0),limits = c(min(fcB), max(fcB)), breaks=Pticks, labels=Pticks1label) +
      xlab(xlabel1) +
      ylab("Percent Failure")

    plotout <- plotout + geom_line(data=df2, aes(Xline,Fline, colour = best_fit), linewidth = 0.9, linetype = "dashed") +
      scale_color_manual(name = best_fit_label, values=col_legend[1:length(databystress)],breaks = line_breaks)

    # Now we can turn off confidence bounds if we set it to anything but NULL (RCS 5/27/2025)
    if(is.null(nobounds) == TRUE){
      plotout <- plotout + geom_ribbon(data=df3, aes(x=Xline2up, ymin=Fline2up, ymax=Fline2down, fill = conf_legend), alpha=0.2) +
        scale_fill_manual(name = paste(c(confid*100,"% Confidence Intervals"),collapse = " "),values=col_legend[1:length(databystress)],breaks = conf_breaks)
    }
  }
  # return(outputpp1[[i]][[2]][1])
  if(is.null(MLE_i) == TRUE){
    # Logging
    stress_set <- rep(0,length(databystress)*length(3:length(databystress[[i]][1,])))
    # return(stress_set)
    lambda_set <- rep(0,length(databystress))
    R2_set <- rep(0,length(databystress))
    for(i in 1:length(databystress)){
      if(length(3:length(databystress[[i]][1,])) == 1){
        stress_set[i] <- c(outputpp1[[i]][[1]])
      } else{
        stress_set[i] <- c(outputpp1[[i]][[1]])[1]
        stress_set[i+length(databystress)] <- c(outputpp1[[i]][[1]])[2]
      }
      lambda_set[i] <- outputpp1[[i]][[2]][1]
      R2_set[i] <- outputpp1[[i]][[3]]
    }
    matset <- c(stress_set,lambda_set,SSEtot,R2_set)

    # Produce some output text that summarizes the results
    cat(c("Exponential Least-Squares estimates\n\n"))
    if(length(3:length(databystress[[i]][1,])) == 1){
      print(matrix(matset, nrow = length(databystress), ncol = 4, byrow = FALSE,dimnames = list(rep("",length(databystress)),c("Stress","\U03BB","SSE","R\U00B2"))))
    }
    if(length(3:length(databystress[[i]][1,])) == 2){
      print(matrix(matset, nrow = length(databystress), ncol = 5, byrow = FALSE,dimnames = list(rep("",length(databystress)),c("Stress1","Stress2","\U03BB","SSE","R\U00B2"))))
    }
    cat("\n")
  } else{
    # Logging
    stress_set <- rep(0,length(databystress)*length(3:length(databystress[[i]][1,])))
    lambda_set <- rep(0,length(databystress))
    lambda_setCIlow <- rep(0,length(databystress))
    lambda_setCIhigh <- rep(0,length(databystress))
    loglik_set <- rep(0,length(databystress))
    AIC_set <- rep(0,length(databystress))
    for(i in 1:length(databystress)){
      if(length(3:length(databystress[[i]][1,])) == 1){
        stress_set[i] <- c(outputpp1[[i]][[1]])
      } else{
        stress_set[i] <- outputpp1[[i]][[1]][1]
        stress_set[i+length(databystress)] <- outputpp1[[i]][[1]][2]
      }
      lambda_set[i] <- outputpp1[[i]][[2]][1]
      lambda_setCIlow[i] <- ttfc_MLE_list[[i]][[3]][[1]][1]
      lambda_setCIhigh[i] <- ttfc_MLE_list[[i]][[3]][[1]][2]
      loglik_set[i] <- ttfc_MLE_list[[i]][[4]]
      AIC_set[i] <- ttfc_MLE_list[[i]][[6]]
    }

    # Produce some output text that summarizes the results
    cat(c("Exponential Maximum Likelihood estimates\n\n"))
    if(length(3:length(databystress[[i]][1,])) == 1){
      matset <- c(stress_set,lambda_set,lambda_setCIlow,lambda_setCIhigh,loglik_set,AIC_set)
      print(matrix(matset, nrow = length(databystress), ncol = 6, byrow = FALSE,dimnames = list(rep("",length(databystress)),c("Stress","\U03BB","\U03BB Lower CI","\U03BB Upper CI","loglik","AIC"))))
    }
    if(length(3:length(databystress[[i]][1,])) == 2){
      matset <- c(stress_set,lambda_set,lambda_setCIlow,lambda_setCIhigh,loglik_set,AIC_set)
      print(matrix(matset, nrow = length(databystress), ncol = 7, byrow = FALSE,dimnames = list(rep("",length(databystress)),c("Stress1","Stress2","\U03BB","\U03BB Lower CI","\U03BB Upper CI","loglik","AIC"))))
    }
    cat("\n")
  }
  return(list(output = outputpp1, summary.nonparametric=xiRFblock_list, prob_plot = plotout))
}
