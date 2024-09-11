# Life-Stress Relationship Plot Generator (LSQ Life-Stress)
# Developed by Dr. Reuel Smith, 2021-2024

lifestress.relationplot.LSQ <- function(data,ls,dist,params,S=NULL,L=NULL,Smin=NULL,Smax=NULL,Suse=NULL,therm=1,confid=0.95,Llab="Characteristic Life - L",Slab="Characteristic Stress - S") {
  # Minimum inputs: Original data and stresses, life-stress model, life distribution, life parameters (alpha, mu, mu_t, etc.)
  # Optional inputs: Use stress, min and max stress, confidence, labels
  # Output: Relationship plot
  # RCS (8/2/2024) - Going to include the distribution overlay on each stress level.  Will need to scale these accordingly based on stress level
  # Base the scale on the minimum distance between stress levels

  # Load plotly library for 3D plotting
  library(plotly)
  library(plyr)
  library(ggplot2)

  # Legend colors
  col_legend <- c("red","blue","darkgreen","violet","lightblue","orange","pink","darkblue","lightgreen","yellow","green")

  # Check first that the data has multiple accelerated stress levels
  if(length(checkstress(data))==1) {
    stop('Need more than one stress level to generate relationship plot.')
  }

  # Then group data by type for naming purposes
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
    fulldatabystress <- databystress
  } else{
    databystress<-checkstress(data)[which(dcount[[1]] == 1 & dcount[[2]] == 1)]
    singledat<-checkstress(data)[which(!dcount[[1]] == 1 & dcount[[2]] == 1)]
    fulldatabystress<-c(singledat,databystress)
  }

  # Reorder the S and L vector inputs to match fulldatabystress
  if(is.null(singledat)==FALSE){
    Snew<-c(S[which(S==c(sort.xircstressdata(data)[[3]][1])):length(S)],S[1:(which(S==c(sort.xircstressdata(data)[[3]][1]))-1)])
    Lnew<-c(L[which(S==c(sort.xircstressdata(data)[[3]][1])):length(S)],L[1:(which(S==c(sort.xircstressdata(data)[[3]][1]))-1)])
    S<-Snew
    L<-Lnew
  }


  # return(c(databystress,singledat))


  # Compute min and max stress used for relationship plot (if single-stress and if not given)
  if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    # TRAIN SMIN FOR CASES WHERE LEFT OUT
    if(is.null(Smin)==TRUE && is.null(Suse)==TRUE){
      if(log10(min(data[,3])) > 1){
        Smin <- round_any(min(data[,3]), 0.5*(10^floor(log10(min(data[,3])))), f = floor)
      }
      # Scale between 1 and 10 for single digit stress
      if(log10(min(data[,3])) < 1 && log10(min(data[,3])) > 0){
        Smin <- 1
      }
      # Scale between .1 and 1 for single digit stress
      if(log10(min(data[,3])) < 0 && log10(min(data[,3])) > -1){
        Smin <- 0.1
      }
    }
    if(is.null(Smin)==TRUE && is.null(Suse)==FALSE && Suse < max(data[,3])){
      if(log10(Suse) > 1){
        Smin <- round_any(Suse, 0.5*(10^floor(log10(Suse))), f = floor)
      }
      # Scale between 1 and 10 for single digit stress
      if(log10(Suse) < 1 && log10(Suse) > 0){
        Smin <- 1
      }
      # Scale between .1 and 1 for single digit stress
      if(log10(Suse) < 0 && log10(Suse) > -1){
        Smin <- 0.1
      }
      # Scale between .1 and 1 for single digit stress
      if(log10(Suse) < 0){
        Smin <- 10^(floor(log10(Suse)))
      }
    }
    # TRAIN SMAX FOR CASES WHERE LEFT OUT
    if(is.null(Smax)==TRUE){
      if(log10(max(data[,3])) > 1){
        Smax <- round_any(max(data[,3]), 0.5*(10^floor(log10(max(data[,3])))), f = ceiling)
      }
      # Scale between 1 and 10 for single digit stress
      if(log10(max(data[,3])) < 1 && log10(max(data[,3])) > 0){
        Smax <- 10
      }
      # Scale between 0.1 and 1 for single digit stress
      if(log10(max(data[,3])) < 0 && log10(max(data[,3])) > -1){
        Smax <- 1
      }
    }
    # Setup or select the Life-Stress lines based on parameter estimates for theta and your
    # minimum and maximum stress values
    Sline<-linspace(Smin,Smax,100)
    if(dist=="Exponential"){
      Lline<-lifestress.select(ls)[[1]](params,Sline)
    }
    if(dist=="3PWeibull"){
      Lline<-lifestress.select(ls)[[1]](params[3:length(params)],Sline)
    }
    if(dist!="3PWeibull" && dist!="Exponential"){
      Lline<-lifestress.select(ls)[[1]](params[2:length(params)],Sline)
    }
  }
  # For dual-stress life-stress models
  if(ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3"){
    if(is.null(Smin)==TRUE && is.null(Suse)==TRUE){
      Smin<-c(0,0)
      Smin[1] <- round_any(min(data[,3]), 0.5*(10^floor(log10(min(data[,3])))), f = floor)
      Smin[2] <- round_any(min(data[,4]), 0.5*(10^floor(log10(min(data[,4])))), f = floor)
    }
    if(is.null(Smin)==TRUE && is.null(Suse)==FALSE){
      Smin<-c(0,0)
      Smin[1] <- round_any(Suse[1], 0.5*(10^floor(log10(Suse[1]))), f = floor)
      Smin[2] <- round_any(Suse[2], 0.5*(10^floor(log10(Suse[2]))), f = floor)
    }
    if(is.null(Smax)==TRUE){
      Smax<-c(0,0)
      Smax[1] <- round_any(max(data[,3]), 0.5*(10^floor(log10(max(data[,3])))), f = ceiling)
      Smax[2] <- round_any(max(data[,4]), 0.5*(10^floor(log10(max(data[,4])))), f = ceiling)
    }
    # Setup or select the Life-Stress function based on parameter estimates for theta and your
    # minimum and maximum stress values
    Sline1<-linspace(Smin[1],Smax[1],100)
    Sline2<-linspace(Smin[2],Smax[2],100)
    if(dist=="Exponential"){
      Lline<-lifestress.select(ls)[[1]](params,Sline1,Sline2)
    }
    if(dist=="3PWeibull"){
      Lline<-lifestress.select(ls)[[1]](params[3:length(params)],Sline1,Sline2)
    }
    if(dist!="3PWeibull" && dist!="Exponential"){
      Lline<-lifestress.select(ls)[[1]](params[2:length(params)],Sline1,Sline2)
    }
  }

  # Form data frame
  if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    # Will likely have to separate these into LSQ, MLE, and Bayesian if I want to have them detect the parameters.
    # Will also need to do this for step-stress and maybe ADT.  Will have this under one help file for simplicity.
    data_legend <- logical(0)
    dist_legend <- logical(0)
    # Lists the data by stress
    for(i in 1:length(fulldatabystress)){
      if(size(fulldatabystress[[i]])[1] > 1){
        data_legend<-c(data_legend,rep(paste(c("Data for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "),sum(fulldatabystress[[i]][,2])))
      } else{
        data_legend<-c(data_legend,paste(c("Data for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
      }
    }
    # Lists the life distribution parameters by stress
    # UPDATE (8/2/2024): Adding distribution bands by stress (maximum and minimum per ggplot)
    dist_scale_S <- min(diff(S))
    dist_scale_invS <- min(diff(1/S))

    for(i in 1:length(fulldatabystress)){
      if(dist == "Normal"){
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
          dist_legend<-c(dist_legend,paste(c("Distribution for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
          dist_legend<-c(dist_legend,paste(c("Distribution for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
        # Pull the life range from the mode to the end
        # xrange <- c()
        # Then compute the lower pdf (to be transformed to ymin)
      }
      if(dist == "Lognormal"){
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
      }
      if(dist == "Exponential"){
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("1/\U03BB for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("1/\U03BB for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
      }
      if(dist == "Weibull"){
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("\U03B1 for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("\U03B1 for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
        # if(params[1] >= 1){
        #   xrangemin <- linspace(L[i]*(((params[1] - 1)/params[1])^(1/params[1])),pweibull(0.999,L[i],params[1]),50)
        #   fmin <- dweibull(xrangemin,L[i],params[1])
        # }
        # if(params[1] < 1){
        #   xrangemin <- linspace(pweibull(0.001,L[i],params[1]),pweibull(0.999,L[i],params[1]),50)
        #   fmin <- dweibull(xrangemin,L[i],params[1])
        # }
      }
      if(dist == "Gumbel"){
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
      }
      if (dist=="Logistic") {
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
      }
      if (dist=="Loglogistic") {
        if(size(fulldatabystress[[i]])[1] > 1){
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][1,3:length(fulldatabystress[[i]][1,])]),collapse = " "))
        } else{
          data_legend<-c(data_legend,paste(c("\U03BC for stress level",fulldatabystress[[i]][3:length(fulldatabystress[[i]])]),collapse = " "))
        }
      }
    }
    if(is.null(Suse) == FALSE){
      if(dist=="Exponential"){
        Luse<-lifestress.select(ls)[[1]](params,Suse)
      }
      if(dist=="3PWeibull"){
        Luse<-lifestress.select(ls)[[1]](params[3:length(params)],Suse)
      }
      if(dist!="3PWeibull" && dist!="Exponential"){
        Luse<-lifestress.select(ls)[[1]](params[2:length(params)],Suse)
      }

      if(dist == "Normal"){
        data_legend<-c(data_legend,paste(c("\U03BC for use stress level",Suse),collapse = " "))
      }
      if(dist == "Lognormal"){
        data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for use stress level",Suse),collapse = " "))
      }
      if(dist == "Exponential"){
        data_legend<-c(data_legend,paste(c("1/\U03BB for use stress level",Suse),collapse = " "))
      }
      if(dist == "Weibull"){
        data_legend<-c(data_legend,paste(c("\U03B1 for use stress level",Suse),collapse = " "))
      }
      if(dist == "Gumbel"){
        data_legend<-c(data_legend,paste(c("\U03BC for use stress level",Suse),collapse = " "))
      }
      if (dist=="Logistic") {
        data_legend<-c(data_legend,paste(c("\U03BC for use stress level",Suse),collapse = " "))
      }
      if (dist=="Loglogistic") {
        data_legend<-c(data_legend,paste(c("\U03BC for use stress level",Suse),collapse = " "))
      }
      df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]],S,Suse), Sinv = 1/c(sort.xircstressdata(data)[[3]],S,Suse), L = c(sort.xircstressdata(data)[[1]],L,Luse),data = data_legend)
    }
    if(is.null(Suse) == TRUE){
      df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]],S), Sinv = 1/c(sort.xircstressdata(data)[[3]],S), L = c(sort.xircstressdata(data)[[1]],L),data = data_legend)
      # df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]],S), Sinv = 1/c(sort.xircstressdata(data)[[3]],S), L = c(sort.xircstressdata(data)[[1]],L))
    }
    df_line <- data.frame(S = Sline, Sinv = 1/Sline, L = Lline, best_fit = rep("Fitted",100))
    # df <- data.frame(X = xrange, YCDF = ycdf, YCDFlow = ycdf_low, YCDFhigh = ycdf_high, best_fit = rep("Fitted",1000))
  }

  # return(list(df_data))

  # UPDATE (11/9/2023): Adding the test of the plots now.  These will be standard output.  Going to test Arrhenius first
  # for all conditions and then do the rest when I'm satisfied.

  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      xlab("Characteristic Stress - S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="Exponential2"){
    # theta[1] - parameter a, theta[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - 1/S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HaS to be in Kelvin for this to work
    K<-8.617385e-5
    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - 1/S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - 1/S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - 1/S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }

  if (ls=="InversePower2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      scale_y_continuous(trans = 'log10') +
      xlab("Characteristic Stress - S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }
  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'blue', size = 0.9, linetype = "dashed") +
      xlab("Characteristic Stress - S") +
      ylab("Characteristic Life - L")

    if(is.null(Suse) == FALSE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)),23)) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
    }
    if(is.null(Suse) == TRUE){
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data, color=data), size=3) +
        scale_shape_manual(values=c(rep(16,length(fulldatabystress)),rep(24,length(fulldatabystress)))) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
    }
  }


  # Check to see that the vectors Smin, Smax, and SUse are the same length to determine
  # if life is based on one stress or multiple stress types
  # if(isTRUE(length(Smin)==length(Smax)&&length(Smin)==length(SUse))){
  #   S_no<-length(Smax)
  # } else {
  #   stop('Make sure that your Smin, Smax and SUse inputs are the same length.')
  # }

  # Sminmaxcheck<-rep(0,S_no)
  #
  # for (i2 in 1:S_no){
  #   if(isTRUE(Smax[i2]>Smin[i2])){
  #     Sminmaxcheck[i2]<-0
  #   } else{
  #     Sminmaxcheck[i2]<-1
  #   }
  # }

  # Check to see if the order of stress max and min is correct (this for one stress)
  # if(sum(Sminmaxcheck)>0) {
  #   stop('Your S_max is less than your S_min.')
  # }



  if(is.null(Suse) == FALSE){
    return(list(df_data,df_line,relationplot=relationplot,Luse=Luse))
    # return(list(relationplot=relationplot,Luse=Luse))
  }
  if(is.null(Suse) == TRUE){
    return(list(df_data,df_line,relationplot=relationplot))
    # return(list(relationplot=relationplot))
  }

  # Sort the stress data for MLE step
  # xircSxiSrc <- sort.xircstressdata(data)

  # # Compute the LSQ and MLE data
  # LSQoutput <- lifestress.LSQest(data,ls,dist,pp)
  # # Check first that the data has multiple accelerated stress levels
  # if(length(LSQoutput[[1]])==1) {
  #   stop('Need more than one stress level to generate relationship plot.')
  # }
  # MLEoutput <- lifestress.MLEest(LSQoutput[[3]],ls,dist,xircSxiSrc[[1]],xircSxiSrc[[3]],xircSxiSrc[[2]],xircSxiSrc[[4]],confid)
  #
  # # ==========================================================================
  # # Here is where we generate the relationship plots
  # # ==========================================================================
  # # Check to see if dist="Exponential" so you can exclude life
  # # distribution parameters.
  # if (dist=="Exponential") {
  #   ishift<-0
  # } else {
  #   ishift<-1
  # }
  #
  # # Setup or select the Life-Stress function based on parameter estimates for theta and your
  # # minimum and maximum stress values
  # if (length(Smax)==1){
  #   Sline<-seq(Smin,Smax,(Smax-Smin)/100)
  # }
  # if (length(Smax)==2){
  #   S1set<-seq(Smin[1],Smax[1],(Smax[1]-Smin[1])/100)
  #   S2set<-seq(Smin[2],Smax[2],(Smax[2]-Smin[2])/100)
  #   S1line <- S1set%o%rep(1,101)
  #   S2line <- rep(1,101)%o%S2set
  #   #S1line <- rep(S1set, each = 100, len = 10100)
  #   #S2line <- rep(S2set,100)
  #   #Sline<-matrix(c(S1line,S2line),nrow=10100,ncol=2,byrow = FALSE)
  # }
  #
  # #
  # # Initialize life-stress parameter estimates for theta
  # if (ls=="Linear") {
  #   # theta[1] - parameter a, theta[2] - parameter b
  #   life <- function(theta,SF) {
  #     theta[ishift+2] + SF*theta[ishift+1]
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Linear Relationship Plot Stress vs. Life
  #   plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=c("Characteristic Stress (",Slab,")"),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
  #   lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(SUse,UselifeLSQ, pch=17, col="blue")
  #   points(SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="Exponential"){
  #   # theta[1] - parameter a, theta[2] - parameter b
  #
  #   life <- function(theta,SF) {
  #     theta[ishift+2]*exp(SF*theta[ishift+1])
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Exponential Relationship Plot Stress vs. Life
  #   plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=c("Characteristic Stress (",Slab,")"),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
  #   lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(SUse,UselifeLSQ, pch=17, col="blue")
  #   points(SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="Arrhenius") {
  #   # lsparams[1] - parameter Ea, lsparams[2] - parameter b
  #   # Temperature HaS to be in Kelvin for this to work
  #   K<-8.617385e-5
  #   life <- function(theta,SF) {
  #     theta[ishift+2]*exp(theta[ishift+1]/(K*SF))
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Arrhenius Relationship Plot Inverse Stress vs. Life
  #   plot(1/Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=expression(paste("Characteristic Stress"^"-1","(","K"^"-1",")")),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(1/Smax,1/Smin),axes=FALSE)
  #   lines(1/Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(1/LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(1/SUse,UselifeLSQ, pch=17, col="blue")
  #   points(1/SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topleft", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="Eyring") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b
  #   life <- function(theta,SF) {
  #     (theta[ishift+2]/SF)*exp(theta[ishift+1]/SF)
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Eyring Relationship Plot Inverse Stress vs. Life
  #   plot(1/Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=expression(paste("Characteristic Stress"^"-1","(","K"^"-1",")")),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(1/Smax,1/Smin),axes=FALSE)
  #   lines(1/Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(1/LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(1/SUse,UselifeLSQ, pch=17, col="blue")
  #   points(1/SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topleft", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="Eyring2") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b
  #   life <- function(theta,SF) {
  #     (1/SF)*exp(-(theta[ishift+1] - (theta[ishift+2]/SF)))
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Eyring 2 Relationship Plot Inverse Stress vs. Life
  #   plot(1/Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=expression(paste("Characteristic Stress"^"-1","(","K"^"-1",")")),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(1/Smax,1/Smin),axes=FALSE)
  #   lines(1/Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(1/LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(1/SUse,UselifeLSQ, pch=17, col="blue")
  #   points(1/SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topleft", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="Power") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b
  #   life <- function(theta,SF) {
  #     theta[ishift+2]*(SF^theta[ishift+1])
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Power Relationship Plot Stress vs. Life
  #   plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=c("Characteristic Stress (",Slab,")"),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
  #   lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(SUse,UselifeLSQ, pch=17, col="blue")
  #   points(SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="InversePower") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b
  #   life <- function(theta,SF) {
  #     theta[ishift+2]*(SF^-theta[ishift+1])
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Inverse Power Relationship Plot Stress vs. Life
  #   plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=c("Characteristic Stress (",Slab,")"),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
  #   lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(SUse,UselifeLSQ, pch=17, col="blue")
  #   points(SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="Logarithmic") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b
  #   life <- function(theta,SF) {
  #     theta[ishift+1]*log(SF) + theta[ishift+2]
  #   }
  #   UselifeLSQ <- life(LSQoutput[[3]],SUse)
  #   UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  #   # Logarithmic Relationship Plot Stress vs. Life
  #   plot(Sline, life(MLEoutput[[1]],Sline), col="black", lty=2,
  #        xlab=c("Characteristic Stress (",Slab,")"),
  #        ylab=c("Characteristic Life (",Llab,")"),
  #        type = "l",pch=13,xlim=c(Smin,Smax),axes=FALSE)
  #   lines(Sline, life(LSQoutput[[3]],Sline), col="blue", lty=3)
  #   points(LSQoutput[[1]],LSQoutput[[2]], pch=16, col="red")
  #   points(SUse,UselifeLSQ, pch=17, col="blue")
  #   points(SUse,UselifeMLE, pch=18, col="black")
  #   axis(1, at=NULL, las=2, cex.axis=0.7)
  #   axis(2, at=NULL, las=2, cex.axis=0.7)
  #   grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = FALSE)
  #   legend(x = "topright", legend=c("Life-Stress Data","LSQ Life-Stress Curve","MLE Life-Stress Curve"),
  #          col=c("red", "blue","black","black"), pch=c(16,-1,-1), lty=c(0,3,2), cex=0.8,
  #          text.font=4, bg='white')
  # }
  #
  # if (ls=="MultiStress") {
  #   # lsparams - c(ao,a1,...,an), S - c(1,S1,...,Sn)
  #   life <- function(theta,SF) {
  #     exp(theta[ishift+1:length(SF)+ishift+1]%*%c(1,SF))
  #   }
  # }
  #
  # if (ls=="TempHumidity") {
  #   # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b
  #   lifeUse <- function(theta,SF) {
  #     theta[ishift+1]*exp((theta[ishift+2]/SF[1]) + (theta[ishift+3]/SF[2]))
  #   } # Single stress entry life (for Use Life primarily)
  #   life <- function(theta,SF1,SF2) {
  #     theta[ishift+1]*exp((theta[ishift+2]/SF1) + (theta[ishift+3]/SF2))
  #   } # Multiple stress entry life (for vectors)
  #   UselifeLSQ <- lifeUse(LSQoutput[[3]],SUse)
  #   UselifeMLE <- lifeUse(MLEoutput[[1]],SUse)
  #
  #   # Temperature-Humidity Relationship Plot Stress vs. Life
  #   persp(x=S1set,y=S2set,z=life(MLEoutput[[1]],S1line,S2line),theta = 100, phi = 50,col="cyan",
  #         shade=0.4, nticks=5, ticktype="detailed",
  #         xlab=c("Characteristic Temperature (",Slab[1],")"),
  #         ylab=c("Relative Humidity (",Slab[2],")"),
  #         zlab=c("Characteristic Life (",Llab,")"))
  #   #fig<-plot_ly(x = ~S1line, y = ~S2line, z = ~life(MLEoutput[[1]],Sline), type = 'mesh3d')
  #   #fig <- fig %>% layout(scene = list(xaxis = list(title = c("Characteristic Stress (",Slab[1],")")),
  #   #                                   yaxis = list(title = c("Characteristic Humidity (",Slab[2],")")),
  #   #                                   zaxis = list(title = c("Characteristic Life (",Llab,")"))))
  #   #fig
  # }
  #
  # if (ls=="TempNonthermal") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
  #   lifeUse <- function(theta,SF) {
  #     theta[ishift+3]/((SF[2]^theta[ishift+2])*exp(-theta[ishift+1]/SF[1]))
  #   } # Single stress entry life (for Use Life primarily)
  #   life <- function(theta,SF1,SF2) {
  #     theta[ishift+3]/((SF2^theta[ishift+2])*exp(-theta[ishift+1]/SF1))
  #   } # Multiple stress entry life (for vectors)
  #   UselifeLSQ <- lifeUse(LSQoutput[[3]],SUse)
  #   UselifeMLE <- lifeUse(MLEoutput[[1]],SUse)
  #
  # }
  #
  # if (ls=="Eyring3") {
  #   # lsparams[1] - parameter a, lsparams[2] - parameter b
  #   # lsparams[3] - parameter c, lsparams[4] - parameter d
  #   life <- function(theta,SF) {
  #     (1/SF[1])*exp((theta[ishift+1] + (theta[ishift+2]/SF[1])) + (theta[ishift+3] + (theta[ishift+4]/SF[1]))*SF[2])
  #   }
  #
  # }
  #
  # # Compute Use life
  # #UselifeLSQ <- life(LSQoutput[[3]],SUse)
  # #UselifeMLE <- life(MLEoutput[[1]],SUse)
  #
  # # Return parameter list
  # return(list(UselifeLSQ,UselifeMLE))
}
