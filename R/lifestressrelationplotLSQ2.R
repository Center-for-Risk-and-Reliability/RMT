# Life-Stress Relationship Plot Generator (LSQ Life-Stress)
# Developed by Dr. Reuel Smith, 2021-2025

lifestress.relationplot.LSQ.2 <- function(data,ls,dist,params,S=NULL,L=NULL,distparams=NULL,
                                        Smin=NULL,Smax=NULL,Suse=NULL,therm=1,confid=0.95,
                                        Llab=NULL,Slab=NULL,Slab2=NULL,confid_int = NULL,predic_int = NULL,
                                        stressunit1 = NULL, stressunit2 = NULL) {
  # Minimum inputs: Original data and stresses, life-stress model, life distribution, life parameters (alpha, mu, mu_t, etc.)
  # Optional inputs: Use stress, min and max stress, confidence, labels
  # Output: Relationship plot
  # UPDATE (11/9/2023): Adding the test of the plots now.  These will be standard output.  Going to test Arrhenius first
  # for all conditions and then do the rest when I'm satisfied.
  # UPDATE RCS (8/2/2024) - Going to include the distribution overlay on each stress level.  Will need to scale these accordingly based on stress level
  # Base the scale on the minimum distance between stress levels
  # UPDATE (8/2/2024): Adding distribution bands by stress (maximum and minimum per ggplot)
  # UPDATE (12/20/2025): (1) Updated the output plots to include stress units like the probability plots.  (2) Reactivated the
  # two-stress life models TempHumidity, TempNonthermal, and Eyring3 to plot dual relationship plots


  # Load plotly library for 3D plotting
  library(plotly)
  library(plyr)
  library(ggplot2)

  # Legend colors
  col_legend <- c("red","blue","darkgreen","violet","gold","orange","pink2","darkblue","lightgreen","yellow","green","darkviolet","darkorange","darkred","purple","royalblue","brown","lightpink","tan","darkgray","aquamarine","sienna","limegreen","mediumpurple3","chocolate","red4")
  # Legend shapes
  shape_legend <- c(0:10)
  # Legend shapes
  shape_legend2 <- c(15:25)
  # Legend line type
  linetype_legend <- rep(c(1,2,4,5,6),5)

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

  # return(list(fulldatabystress,singledat))
  # Reorder the S and L vector inputs to match fulldatabystress
  if(is.null(singledat)==FALSE){ # sort S, L, and distparams differently if there is single data included
    S.new <- S
    # L.new <- rep(0,length(S))
    S.new<-S
    L.new<-L
    # S.new<-c(S[which(S==c(sort.xircstressdata(data)[[3]][1])):length(S)],S[1:(which(S==c(sort.xircstressdata(data)[[3]][1]))-1)])
    # L.new<-c(L[which(S==c(sort.xircstressdata(data)[[3]][1])):length(S)],L[1:(which(S==c(sort.xircstressdata(data)[[3]][1]))-1)])
    # distparams.new <- c()
    S<-S.new
    L<-L.new
  }
  # return(list(S,L,distparams,fulldatabystress))

  # return(list(sort(S),S[which(S==c(sort.xircstressdata(data)[[3]][1])):length(S)],S[1:(which(S==c(sort.xircstressdata(data)[[3]][1]))-1)]))
  # Form data frame
  # Will likely have to separate these into LSQ, MLE, and Bayesian if I want to have them detect the parameters.
  # Will also need to do this for step-stress and maybe ADT.  Will have this under one help file for simplicity.
  data_legend <- logical(0)
  dist_legend <- logical(0)
  data_breaks <- rep(0,length(fulldatabystress))
  pdf_breaks <- rep(0,length(fulldatabystress))
  # return(fulldatabystress)
  # Lists the data by stress
  for(i in 1:length(fulldatabystress)){
    # UPDATED FORMAT FOR DATA LEGEND CREATION
    if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
      # Only plot data if multiple
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        if(length(databystress) > 1){ # multiple groups of data
          data_legend<-c(data_legend,rep(paste(c("Data for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),sum(fulldatabystress[[i]][,2])))
          data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
        } else{ # only one group
          data_legend<-c(data_legend,rep("Data",sum(fulldatabystress[[i]][,2])))
          data_breaks <- waiver()
        }
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,rep(paste(c("Data for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),sum(fulldatabystress[[i]][,2])))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),sum(fulldatabystress[[i]][,2])))
        data_breaks[i] <- paste(c("Data for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),sum(fulldatabystress[[i]][,2])))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),sum(fulldatabystress[[i]][,2])))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,rep(paste(c("Data for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),sum(fulldatabystress[[i]][,2])))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
      }
    }
    # *****************************************
    if(is.vector(fulldatabystress[[i]])==TRUE){ # Indicate that the block is a single-data vector
      if(length(3:length(fulldatabystress[[i]])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        if(length(databystress) > 1){ # multiple groups of data
          data_legend<-c(data_legend,paste(c("Data for",fulldatabystress[[i]][3]," units"),collapse = " "))
          data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][3]," units"),collapse = " ")
        } else{ # only one group
          data_legend<-c(data_legend,"Data")
          data_breaks <- waiver()
        }
      }
      if(length(3:length(fulldatabystress[[i]])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("Data for",fulldatabystress[[i]][3]," ",stressunit1),collapse = " "))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][3]," ",stressunit1),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("Data for ",fulldatabystress[[i]][3]," units/",fulldatabystress[[i]][4]," units"),collapse = " "))
        data_breaks[i] <- paste(c("Data for ",fulldatabystress[[i]][3]," units/",fulldatabystress[[i]][4]," units"),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("Data for",fulldatabystress[[i]][3],"",stressunit1,"/",fulldatabystress[[i]][4],"",stressunit2),collapse = " "))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][3],"",stressunit1,"/",fulldatabystress[[i]][4],"",stressunit2),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("Data for",fulldatabystress[[i]][3],"",stressunit1,"/",fulldatabystress[[i]][4]," units"),collapse = " "))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][3],"",stressunit1,"/",fulldatabystress[[i]][4]," units"),collapse = " ")
      }
      if(length(3:length(fulldatabystress[[i]])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("Data for",fulldatabystress[[i]][3],""," units","/",fulldatabystress[[i]][4],"",stressunit2),collapse = " "))
        data_breaks[i] <- paste(c("Data for",fulldatabystress[[i]][3]," units","/",fulldatabystress[[i]][4],"",stressunit2),collapse = " ")
      }
    }
  }
  # return(list(data_legend,L,distparams,fulldatabystress))

  # Lists the life distribution parameters by stress
  # UPDATE (8/2/2024): Adding distribution bands by stress (maximum and minimum per ggplot)
  if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" ||
      ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    dist_scale_S <- abs(min(diff(S)))
    dist_scale_invS <- abs(min(diff(1/S)))
  }
  if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
    dist_scale_S1 <- min(abs(diff(S[,1]))[which(abs(diff(S[,1]))>0)])
    dist_scale_invS1 <- min(abs(diff(1/S[,1]))[which(abs(diff(1/S[,1]))>0)])
    dist_scale_S2 <- min(abs(diff(S[,2]))[which(abs(diff(S[,2]))>0)])
    dist_scale_invS2 <- min(abs(diff(1/S[,2]))[which(abs(diff(1/S[,2]))>0)])
  }

  # Initialize dist.plot.x and dist.plot.y as zero length objects
  dist.plot.x <- numeric(0)
  dist.plot.y <- numeric(0)
  L.cut <- numeric(0)
  S.cut <- numeric(0)
  # ============================================================================
  # Log and sort all of the data points by stress level.  The distribution based
  # life is also sorted but only for cases of multiple data per stress level.
  # ============================================================================
  for(i in 1:length(fulldatabystress)){
    # As this sorts through fulldatabystress, it needs to check and see that it is not a vector (#1),
    # then it needs to match the stress fulldatabystress[[i]][3] to S to pull the distribution parameter (#2)
    if(dist == "Normal"){
      # Generate band for All Stresses
      x.band <- linspace(qnorm(0.01,L[i],distparams[i]),qnorm(0.99,L[i],distparams[i]),100)        # Life Axis
      y.band <- dnorm(x.band,L[i],distparams[i])                                                   # Stress Axis
      if(i == 1){ # Builds full distribution generation table and associated stresses
        dist.plot.x <- x.band
        dist.plot.y <- y.band
      } else{
        dist.plot.x <- c(dist.plot.x,x.band)
        dist.plot.y <- c(dist.plot.y,y.band)
      }

      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
          if(length(databystress) > 1){ # multiple groups of data
            data_legend <- c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
          } else{ # only one group
            data_legend<-c(data_legend,"Data")
            dist_legend <- c(dist_legend,"PDF")
            pdf_breaks[i] <- waiver()
          }
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
          data_legend <- c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
        }
        # DUAL STRESS OPTIONS TO FOLLOW
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
          data_legend <- c(data_legend,paste(c("\U03BC for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
      }
    }

    if(dist == "Lognormal"){
      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix (RCS 1/6/2026)
        # Pull distribution parameter from distparams
        distparams.pull <- distparams[which(S==fulldatabystress[[i]][1,3])]
        # Generate band for All Stresses
        x.band <- linspace(qlnorm(0.01,log(L[which(S==fulldatabystress[[i]][1,3])]),distparams.pull),qlnorm(0.99,log(L[which(S==fulldatabystress[[i]][1,3])]),distparams.pull),100)  # Life Axis
        y.band <- dlnorm(x.band,log(L[which(S==fulldatabystress[[i]][1,3])]),distparams.pull)                                                   # Stress Axis
        # return(list(distparams.pull,x.band,y.band,i))
        if(length(x.band) == 0){ # Builds full distribution generation table and associated stresses
          dist.plot.x <- x.band
          dist.plot.y <- y.band
          L.cut <- L[which(S==fulldatabystress[[i]][1,3])]
          S.cut <- S[which(S==fulldatabystress[[i]][1,3])]
        } else{
          dist.plot.x <- c(dist.plot.x,x.band)
          dist.plot.y <- c(dist.plot.y,y.band)
          L.cut <- c(L.cut,L[which(S==fulldatabystress[[i]][1,3])])
          S.cut <- c(S.cut,S[which(S==fulldatabystress[[i]][1,3])])
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
          if(length(databystress) > 1){ # multiple groups of data
            data_legend <- c(data_legend,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
          } else{ # only one group
            data_legend<-c(data_legend,"Data")
            dist_legend <- c(dist_legend,"PDF")
            pdf_breaks[i] <- waiver()
          }
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
          data_legend <- c(data_legend,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
        }
        # DUAL STRESS OPTIONS TO FOLLOW
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
          data_legend <- c(data_legend,paste(c("exp(\U03BC_t) for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
          data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
          data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
          data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
      }
    }
    if(dist == "Exponential"){
      # Generate band for All Stresses
      x.band <- linspace(qexp(0.01,1/L[i]),qexp(0.99,1/L[i]),100)        # Life Axis
      y.band <- dexp(x.band,1/L[i])                                                   # Stress Axis
      if(i == 1){ # Builds full distribution generation table and associated stresses
        dist.plot.x <- x.band
        dist.plot.y <- y.band
      } else{
        dist.plot.x <- c(dist.plot.x,x.band)
        dist.plot.y <- c(dist.plot.y,y.band)
      }

      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
          if(length(databystress) > 1){ # multiple groups of data
            data_legend <- c(data_legend,paste(c("1/\U03BB for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
          } else{ # only one group
            data_legend<-c(data_legend,"Data")
            dist_legend <- c(dist_legend,"PDF")
            pdf_breaks[i] <- waiver()
          }
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
          data_legend <- c(data_legend,paste(c("1/\U03BB for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("1/\U03BB for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
        }
        # DUAL STRESS OPTIONS TO FOLLOW
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
          data_legend <- c(data_legend,paste(c("1/\U03BB for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("1/\U03BB for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
          data_legend<-c(data_legend,paste(c("1/\U03BB for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("1/\U03BB for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
          data_legend<-c(data_legend,paste(c("1/\U03BB for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("1/\U03BB for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
          data_legend<-c(data_legend,paste(c("1/\U03BB for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("1/\U03BB for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
      }
    }
    if(dist == "Weibull"){
      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix (RCS 1/6/2026)
        # Pull distribution parameter from distparams
        distparams.pull <- distparams[which(S==fulldatabystress[[i]][1,3])]
        # Generate band for All Stresses
        x.band <- linspace(qweibull(0.01,distparams.pull,L[which(S==fulldatabystress[[i]][1,3])]),qweibull(0.99,distparams.pull,L[which(S==fulldatabystress[[i]][1,3])]),100)  # Life Axis
        y.band <- dweibull(x.band,distparams.pull,L[which(S==fulldatabystress[[i]][1,3])])                                                  # Stress Axis
        if(length(x.band) == 0){ # Builds full distribution generation table and associated stresses
          dist.plot.x <- x.band
          dist.plot.y <- y.band
          L.cut <- L[which(S==fulldatabystress[[i]][1,3])]
          S.cut <- S[which(S==fulldatabystress[[i]][1,3])]
        } else{
          dist.plot.x <- c(dist.plot.x,x.band)
          dist.plot.y <- c(dist.plot.y,y.band)
          L.cut <- c(L.cut,L[which(S==fulldatabystress[[i]][1,3])])
          S.cut <- c(S.cut,S[which(S==fulldatabystress[[i]][1,3])])
        }

        if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
          if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
            if(length(databystress) > 1){ # multiple groups of data
              data_legend <- c(data_legend,paste(c("\U03B1 for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
              dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
              pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
            } else{ # only one group
              data_legend<-c(data_legend,"Data")
              dist_legend <- c(dist_legend,"PDF")
              pdf_breaks[i] <- waiver()
            }
          }
          if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
            data_legend <- c(data_legend,paste(c("\U03B1 for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
            data_breaks <- c(data_breaks,paste(c("\U03B1 for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
          }
          # DUAL STRESS OPTIONS TO FOLLOW
          if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
            data_legend <- c(data_legend,paste(c("\U03B1 for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
            data_breaks <- c(data_breaks,paste(c("\U03B1 for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
            pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
          }
          if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
            data_legend<-c(data_legend,paste(c("\U03B1 for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
            data_breaks <- c(data_breaks,paste(c("\U03B1 for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
          }
          if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
            data_legend<-c(data_legend,paste(c("\U03B1 for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
            data_breaks <- c(data_breaks,paste(c("\U03B1 for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
          }
          if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
            data_legend<-c(data_legend,paste(c("\U03B1 for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
            data_breaks <- c(data_breaks,paste(c("\U03B1 for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
          }
        }
      }
    }
    if(dist == "Gumbel"){
      # Generate band for All Stresses
      x.band <- linspace(L[i]+distparams[i]*log(-log(1-0.01)),L[i]+distparams[i]*log(-log(1-0.99)),100)        # Life Axis
      y.band <- (1/distparams[i])*exp(((x.band - L[i])/distparams[i]) - exp((x.band - L[i])/distparams[i]))                                                   # Stress Axis
      if(i == 1){ # Builds full distribution generation table and associated stresses
        dist.plot.x <- x.band
        dist.plot.y <- y.band
      } else{
        dist.plot.x <- c(dist.plot.x,x.band)
        dist.plot.y <- c(dist.plot.y,y.band)
      }

      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
          if(length(databystress) > 1){ # multiple groups of data
            data_legend <- c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
          } else{ # only one group
            data_legend<-c(data_legend,"Data")
            dist_legend <- c(dist_legend,"PDF")
            pdf_breaks[i] <- waiver()
          }
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
          data_legend <- c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
        }
        # DUAL STRESS OPTIONS TO FOLLOW
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
          data_legend <- c(data_legend,paste(c("\U03BC for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
      }
    }
    if (dist=="Logistic") {
      # Generate band for All Stresses
      x.band <- linspace(qlogis(0.01,L[i],distparams[i]),qlogis(0.99,L[i],distparams[i]),100)        # Life Axis
      y.band <- dlogis(x.band,L[i],distparams[i])                                                   # Stress Axis
      if(i == 1){ # Builds full distribution generation table and associated stresses
        dist.plot.x <- x.band
        dist.plot.y <- y.band
      } else{
        dist.plot.x <- c(dist.plot.x,x.band)
        dist.plot.y <- c(dist.plot.y,y.band)
      }

      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
          if(length(databystress) > 1){ # multiple groups of data
            data_legend <- c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
          } else{ # only one group
            data_legend<-c(data_legend,"Data")
            dist_legend <- c(dist_legend,"PDF")
            pdf_breaks[i] <- waiver()
          }
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
          data_legend <- c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
        }
        # DUAL STRESS OPTIONS TO FOLLOW
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
          data_legend <- c(data_legend,paste(c("\U03BC for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
          data_legend<-c(data_legend,paste(c("\U03BC for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("\U03BC for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
      }
    }
    if (dist=="Loglogistic") {
      # Generate band for All Stresses
      x.band <- linspace(exp(log(L[i]) + distparams[i]*log((1/0.99) - 1)),exp(log(L[i]) + distparams[i]*log((1/0.01) - 1)),100)  # Life Axis
      y.band <- exp(((log(x.band) - log(L[i]))/distparams[i]))/(distparams[i]*x.band*(1+exp(((log(x.band) - log(L[i]))/distparams[i])))^2)
      if(i == 1){ # Builds full distribution generation table and associated stresses
        dist.plot.x <- x.band
        dist.plot.y <- y.band
      } else{
        dist.plot.x <- c(dist.plot.x,x.band)
        dist.plot.y <- c(dist.plot.y,y.band)
      }

      if(is.vector(fulldatabystress[[i]])==FALSE){ # Indicate that the block is a multi-data matrix
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
          if(length(databystress) > 1){ # multiple groups of data
            data_legend <- c(data_legend,paste(c("exp(\U03BC) for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "))
            dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " "),100))
            pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," units"),collapse = " ")
          } else{ # only one group
            data_legend<-c(data_legend,"Data")
            dist_legend <- c(dist_legend,"PDF")
            pdf_breaks[i] <- waiver()
          }
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
          data_legend <- c(data_legend,paste(c("exp(\U03BC) for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][[1,3:length(fulldatabystress[[i]][1,])]]," ",stressunit1),collapse = " ")
        }
        # DUAL STRESS OPTIONS TO FOLLOW
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
          data_legend <- c(data_legend,paste(c("exp(\U03BC) for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for ",fulldatabystress[[i]][1,3]," units/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
          data_legend<-c(data_legend,paste(c("exp(\U03BC) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
          data_legend<-c(data_legend,paste(c("exp(\U03BC) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3],"",stressunit1,"/",fulldatabystress[[i]][1,4]," units"),collapse = " ")
        }
        if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
          data_legend<-c(data_legend,paste(c("exp(\U03BC) for",fulldatabystress[[i]][1,3],""," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          dist_legend <- c(dist_legend,rep(paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "),100))
          data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " "))
          pdf_breaks[i] <- paste(c("PDF for",fulldatabystress[[i]][1,3]," units","/",fulldatabystress[[i]][1,4],"",stressunit2),collapse = " ")
        }
      }
    }
  }
  # ============================================================================
  # Set up stress shift after defining the data and distribution for the plot (1/6/2026)
  # ============================================================================
  for(i in 1:length(fulldatabystress)){ # Build the stress shift according to the life-stress model
    if(i == 1){
      if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" ||
          ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
        if(is.null(singledat)==TRUE){ # Default all multiple data
          x.shift.S <- rep(S[i],100)
          x.shift.Sinv <- rep(1/S[i],100)
        }
        if(is.null(singledat)==FALSE && i <= length(S.cut)){ # For single data
          x.shift.S <- rep(S.cut[i],100)
          x.shift.Sinv <- rep(1/S.cut[i],100)
        }
      }
      if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
        if(is.null(singledat)==TRUE){ # Default all multiple data
          x.shift.S1 <- rep(S[i,1],100)
          x.shift.Sinv1 <- rep(1/S[i,1],100)
          x.shift.S2 <- rep(S[i,2],100)
          x.shift.Sinv2 <- rep(1/S[i,2],100)
        }
        if(is.null(singledat)==FALSE  && i <= length(S.cut)){ # For single data (NOTE RCS1/6/2026: Full implementing of dual stress for single data cases need to be finalized)
          x.shift.S1 <- rep(S.cut[i,1],100)
          x.shift.Sinv1 <- rep(1/S.cut[i,1],100)
          x.shift.S2 <- rep(S.cut[i,2],100)
          x.shift.Sinv2 <- rep(1/S.cut[i,2],100)
        }
      }
    } else{
      if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" ||
          ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){

        if(is.null(singledat)==TRUE){ # Default all multiple data
          x.shift.S <- c(x.shift.S,rep(S[i],100))
          x.shift.Sinv <- c(x.shift.Sinv,rep(1/S[i],100))
        }
        if(is.null(singledat)==FALSE  && i <= length(S.cut)){ # For single data
          x.shift.S <- c(x.shift.S,rep(S.cut[i],100))
          x.shift.Sinv <- c(x.shift.Sinv,rep(1/S.cut[i],100))
        }
      }
      if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
        if(is.null(singledat)==TRUE){ # Default all multiple data
          x.shift.S1 <- c(x.shift.S1,rep(S[i,1],100))
          x.shift.Sinv1 <- c(x.shift.Sinv1,rep(1/S[i,1],100))
          x.shift.S2 <- c(x.shift.S2,rep(S[i,2],100))
          x.shift.Sinv2 <- c(x.shift.Sinv2,rep(1/S[i,2],100))
        }
        if(is.null(singledat)==FALSE){ # For single data
          x.shift.S1 <- c(x.shift.S1,rep(S.cut[i,1],100))
          x.shift.Sinv1 <- c(x.shift.Sinv1,rep(1/S.cut[i,1],100))
          x.shift.S2 <- c(x.shift.S2,rep(S.cut[i,2],100))
          x.shift.Sinv2 <- c(x.shift.Sinv2,rep(1/S.cut[i,2],100))
        }
      }
    }
  }
  # ============================================================================
  # To include the use stress part of the plot
  # ============================================================================
  if(is.null(Suse) == FALSE){
    if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" ||
        ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
      if(dist=="Exponential"){
        Luse<-lifestress.select(ls)[[1]](params,Suse)
      }
      if(dist=="3PWeibull"){
        Luse<-lifestress.select(ls)[[1]](params[3:length(params)],Suse)
      }
      if(dist!="3PWeibull" && dist!="Exponential"){
        Luse<-lifestress.select(ls)[[1]](params[2:length(params)],Suse)
      }

      x.shift.S <- c(x.shift.S,rep(Suse,100))
      x.shift.Sinv <- c(x.shift.Sinv,rep(1/Suse,100))
    }
    if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
      if(dist=="Exponential"){
        Luse<-lifestress.select(ls)[[1]](params,Suse[1],Suse[2])
      }
      if(dist=="3PWeibull"){
        Luse<-lifestress.select(ls)[[1]](params[3:length(params)],Suse[1],Suse[2])
      }
      if(dist!="3PWeibull" && dist!="Exponential"){
        Luse<-lifestress.select(ls)[[1]](params[2:length(params)],Suse[1],Suse[2])
      }
      x.shift.S1 <- c(x.shift.S1,rep(Suse[1],100))
      x.shift.Sinv1 <- c(x.shift.Sinv1,rep(1/Suse[1],100))
      x.shift.S2 <- c(x.shift.S2,rep(Suse[2],100))
      x.shift.Sinv2 <- c(x.shift.Sinv2,rep(1/Suse[2],100))
    }

    if(dist == "Normal"){
      # Generate band for Use Stress
      x.use <- linspace(qnorm(0.01,Luse,params[1]),qnorm(0.99,Luse,params[1]),100) # Life Axis
      y.use <- dnorm(x.use,Luse,params[1])                                               # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("\U03BC for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("\U03BC for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend <- c(data_legend,paste(c("\U03BC for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("\U03BC for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }
    if(dist == "Lognormal"){
      # Generate band for Use Stress
      x.use <- linspace(qlnorm(0.01,log(Luse),params[1]),qlnorm(0.99,log(Luse),params[1]),100) # Life Axis
      y.use <- dlnorm(x.use,log(Luse),params[1])                                               # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("exp(\U03BC_t) for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("exp(\U03BC_t) for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("exp(\U03BC_t) for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("exp(\U03BC_t) for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC_t) for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }
    if(dist == "Exponential"){
      # Generate band for Use Stress
      x.use <- linspace(qexp(0.01,1/Luse),qexp(0.99,1/Luse),100) # Life Axis
      y.use <- dexp(x.use,1/Luse)                                               # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("1/\U03BB for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("1/\U03BB for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("1/\U03BB for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("1/\U03BB for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("1/\U03BB for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("1/\U03BB for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("1/\U03BB for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("1/\U03BBC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("1/\U03BB for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("1/\U03BB for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("1/\U03BB for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("1/\U03BB for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }
    if(dist == "Weibull"){
      # Generate band for Use Stress
      x.use <- linspace(qweibull(0.01,params[1],Luse),qweibull(0.99,params[1],Luse),100) # Life Axis
      y.use <- dweibull(x.use,params[1],Luse)                                               # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("\U03B1 for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("\U03B1 for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("\U03B1 for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03B1 for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("\U03B1 for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03B1 for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03B1 for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03B1 for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03B1 for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03B1 for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03B1 for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03B1 for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }
    if(dist == "Gumbel"){
      # Generate band for Use Stress
      x.use <- linspace(Luse+params[1]*log(-log(1-0.01)),Luse+params[1]*log(-log(1-0.99)),100)        # Life Axis
      y.use <- (1/params[1])*exp(((x.use - Luse)/params[1]) - exp((x.use - Luse)/params[1]))          # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("\U03BC for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("\U03BC for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("\U03BC for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("\U03BC for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }
    if (dist=="Logistic") {
      # Generate band for Use Stress
      x.use <- linspace(qlogis(0.01,Luse,params[1]),qlogis(0.99,Luse,params[1]),100)                                           # Life Axis
      y.use <- dlogis(x.use,Luse,params[1])                                                                                    # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("\U03BC for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("\U03BC for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("\U03BC for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("\U03BC for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("\U03BC for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("\U03BC for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }
    if (dist=="Loglogistic") {
      # Generate band for Use Stress
      x.use <- linspace(exp(log(Luse) + params[1]*log((1/0.99) - 1)),exp(log(Luse) + params[1]*log((1/0.01) - 1)),100)         # Life Axis
      y.use <- exp(((log(x.use) - log(Luse))/params[1]))/(params[1]*x.use*(1+exp(((log(x.use) - log(Luse))/params[1])))^2)     # Stress Axis
      dist.plot.x <- c(dist.plot.x,x.use)
      dist.plot.y <- c(dist.plot.y,y.use)

      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == TRUE){ # When no units then we can negate stress designation
        data_legend <- c(data_legend,paste(c("exp(\U03BC) for use stress",Suse," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," units"),collapse = " "),100))
        data_breaks <- c(pdf_breaks,paste(c("exp(\U03BC) for use stress",Suse," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 1 && is.null(stressunit1) == FALSE && (is.null(stressunit2) == TRUE || is.null(stressunit2) == FALSE)){ # When units are given, we state stress and units even when one group is given
        data_legend<-c(data_legend,paste(c("exp(\U03BC) for use stress",Suse," ",stressunit1),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for use stress",Suse," ",stressunit1),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for use stress",Suse," ",stressunit1),collapse = " "))
      }
      # DUAL STRESS OPTIONS TO FOLLOW
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == TRUE){  # When no units are given and we have two stresses
        data_legend <- c(data_legend,paste(c("exp(\U03BC) for use stress ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for ",Suse[1]," units/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == FALSE){ # When units are given and we have two stresses
        data_legend<-c(data_legend,paste(c("exp(\U03BC) for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2],"",stressunit2),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == FALSE && is.null(stressunit2) == TRUE){ # When units are given for ONLY stress 1 and we have two stresses
        data_legend<-c(data_legend,paste(c("exp(\U03BC) for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1],"",stressunit1,"/",Suse[2]," units"),collapse = " "))
      }
      if(length(3:length(fulldatabystress[[i]][1,])) == 2 && is.null(stressunit1) == TRUE && is.null(stressunit2) == FALSE){ # When units are given for ONLY stress 2 and we have two stresses
        data_legend<-c(data_legend,paste(c("exp(\U03BC) for",Suse[1],""," units","/",Suse[2],"",stressunit2),collapse = " "))
        dist_legend <- c(dist_legend,rep(paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "),100))
        data_breaks <- c(data_breaks,paste(c("exp(\U03BC) for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
        pdf_breaks <- c(pdf_breaks,paste(c("PDF for",Suse[1]," units","/",Suse[2],"",stressunit2),collapse = " "))
      }
    }

    # return(list(sort.xircstressdata(data)[[3]],S.cut,Suse,sort.xircstressdata(data)[[1]],L.cut,Luse,data_legend))
    # SET UP THE DATA IN THE CASE OF HAVING SUSE DEFINED
    if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" ||
        ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
      if(is.null(singledat)==TRUE){ # Default all multiple data
        df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]],S,Suse), Sinv = 1/c(sort.xircstressdata(data)[[3]],S,Suse), L = c(sort.xircstressdata(data)[[1]],L,Luse),data.group = data_legend)
      }
      if(is.null(singledat)==FALSE){ # For single data
        df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]],S.cut,Suse), Sinv = 1/c(sort.xircstressdata(data)[[3]],S.cut,Suse), L = c(sort.xircstressdata(data)[[1]],L.cut,Luse),data.group = data_legend)
      }
    }
    if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
      if (ls=="TempHumidity") {
        if(is.null(Slab)==TRUE){
          Slab <- "Characteristic Temperature - S"
        }
        if(is.null(Slab2)==TRUE){
          Slab2 <- "Characteristic Humidity - H"
        }
        if(is.null(Llab)==TRUE){
          Llab <- "Characteristic Life - L"
        }
      }
      if (ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4") {
        if(is.null(Slab)==TRUE){
          Slab <- "Characteristic Temperature - S"
        }
        if(is.null(Slab2)==TRUE){
          Slab2 <- "Characteristic Nonthermal Stress - U"
        }
        if(is.null(Llab)==TRUE){
          Llab <- "Characteristic Life - L"
        }
      }
      df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]][,1],S[,1],Suse[1],sort.xircstressdata(data)[[3]][,2],S[,2],Suse[2]),
                            Sinv = c(sort.xircstressdata(data)[[3]][,1],S[,1],Suse[1],sort.xircstressdata(data)[[3]][,2],S[,2],Suse[2]),
                            L = rep(c(sort.xircstressdata(data)[[1]],L,Luse)),
                            data.group = rep(data_legend,2),
                            stress.group = c(rep(Slab,length(c(sort.xircstressdata(data)[[3]][,1],S[,1],Suse[1]))),rep(Slab2,length(c(sort.xircstressdata(data)[[3]][,1],S[,1],Suse[1])))))
    }
  }
  # return(list(sort.xircstressdata(data)[[3]],S,Suse,sort.xircstressdata(data)[[1]],L,Luse,df_data))

  if(is.null(Suse) == TRUE){
    # SET UP THE DATA IN THE CASE OF HAVING NO SUSE DEFINED
    if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" ||
        ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
      df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]],S), Sinv = 1/c(sort.xircstressdata(data)[[3]],S), L = c(sort.xircstressdata(data)[[1]],L),data.group = data_legend)
    }
    if (ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
      if (ls=="TempHumidity") {
        if(is.null(Slab)==TRUE){
          Slab <- "Characteristic Temperature - S"
        }
        if(is.null(Slab2)==TRUE){
          Slab2 <- "Characteristic Humidity - H"
        }
        if(is.null(Llab)==TRUE){
          Llab <- "Characteristic Life - L"
        }
      }
      if (ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4") {
        if(is.null(Slab)==TRUE){
          Slab <- "Characteristic Temperature - S"
        }
        if(is.null(Slab2)==TRUE){
          Slab2 <- "Characteristic Nonthermal Stress - U"
        }
        if(is.null(Llab)==TRUE){
          Llab <- "Characteristic Life - L"
        }
      }

      df_data <- data.frame(S = c(sort.xircstressdata(data)[[3]][,1],S[,1],sort.xircstressdata(data)[[3]][,2],S[,2]),
                            Sinv = c(sort.xircstressdata(data)[[3]][,1],S[,1],sort.xircstressdata(data)[[3]][,2],S[,2]),
                            L = rep(c(sort.xircstressdata(data)[[1]],L)),
                            data.group = rep(data_legend,2),
                            stress.group = c(rep(Slab,length(c(sort.xircstressdata(data)[[3]][,1],S[,1]))),rep(Slab2,length(c(sort.xircstressdata(data)[[3]][,1],S[,1])))))
    }
  }

  # Best fit line Analysis
  # Compute min and max stress used for relationship plot (if single-stress and if not given)
  if (ls=="Linear" || ls=="Exponential" || ls=="Exponential2" || ls=="Arrhenius" || ls=="Eyring" || ls=="Eyring2" || ls=="Power" || ls=="PowerwithBias" || ls=="InversePower" || ls=="InversePower2" || ls=="InversePower2" || ls=="Logarithmic"){
    # TRAIN SMIN FOR CASES WHERE LEFT OUT
    if(is.null(Smin)==TRUE && is.null(Suse)==TRUE){
      Smin <- 0
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
      Smin <- 0
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
    # if(is.null(Suse)==FALSE && (ls == "Exponential2")){
    #   Smin <- 1/((1/Suse) + dist_scale_invS*0.8)
    # }
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
    df_line <- data.frame(S = Sline, Sinv = 1/Sline, L = Lline, best_fit = rep("Fitted",100))
  }
  # For dual-stress life-stress models
  # This will need specific life-stress lines as one of the stresses need to be held constant
  # Also Humidity needs to be held between 0.001 and 1
  if(ls=="TempHumidity" || ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
    if(is.null(Smin)==TRUE && is.null(Suse)==TRUE){ # No Smin or Suse defined
      Smin<-c(0,0)
      Smin[1] <- round_any(min(data[,3]), 0.5*(10^floor(log10(min(data[,3])))), f = floor)
      if (ls=="TempHumidity"){
        Smin[2] <- round_any(min(data[,4]), 0.5*(10^floor(log10(min(data[,4])))), f = floor)-0.1
      } else{
        Smin[2] <- round_any(min(data[,4]), 0.5*(10^floor(log10(min(data[,4])))), f = floor)
      }
    }
    if(is.null(Smin)==TRUE && is.null(Suse)==FALSE){ # No Smin defined but there is an Suse
      Smin<-c(0,0)
      Smin[1] <- round_any(Suse[1], 0.5*(10^floor(log10(Suse[1]))), f = floor)
      if (ls=="TempHumidity"){
        Smin[2] <- round_any(Suse[2], 0.5*(10^floor(log10(Suse[2]))), f = floor)-0.1
      } else{
        Smin[2] <- round_any(Suse[2], 0.5*(10^floor(log10(Suse[2]))), f = floor)
      }
    }
    if(is.null(Smax)==TRUE){ # no defined Smax
      Smax<-c(0,0)
      Smax[1] <- round_any(max(data[,3]), 0.5*(10^floor(log10(max(data[,3])))), f = ceiling)
      if (ls=="TempHumidity"){
        Smax[2]<-1
      } else{
        Smax[2] <- round_any(max(data[,4]), 0.5*(10^floor(log10(max(data[,4])))), f = ceiling)
      }
    }
    # Setup or select the Life-Stress function based on parameter estimates for theta and your
    # minimum and maximum stress values
    Sline1.0<-linspace(Smin[1],Smax[1],100)
    Sline2.0<-linspace(Smin[2],Smax[2],100)
    if(ls=="TempHumidity"){
      # First plot line with changing temperature
      if(is.null(Suse)==FALSE){
        for(i in 1:length(unique(rbind(S,Suse)[,1]))){
          if(dist=="Exponential"){
            # unique(rbind(S,Suse)[,1])[i]
            # unique(rbind(S,Suse)[,2])[i]
            Lline1.0<-lifestress.select(ls)[[1]](params,Sline1.0,unique(rbind(S,Suse)[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params,unique(rbind(S,Suse)[,1])[i],Sline2.0)
          }
          if(dist!="Exponential"){
            Lline1.0<-lifestress.select(ls)[[1]](params[2:length(params)],Sline1.0,unique(rbind(S,Suse)[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params[2:length(params)],unique(rbind(S,Suse)[,1])[i],Sline2.0)
          }
          if(i==1){
            Lline1 <- Lline1.0
            Lline2 <- Lline2.0
            Sline1 <- Sline1.0
            Sline2 <- Sline2.0
            Life_legend1 <- rep(paste(c("RH at",unique(rbind(S,Suse)[,2])[i]," ",stressunit2),collapse = " "),100)
            Life_legend2 <- rep(paste(c("Temperature at",unique(rbind(S,Suse)[,1])[i]," ",stressunit1),collapse = " "),100)
          } else{
            Lline1 <- c(Lline1,NA,Lline1.0)
            Lline2 <- c(Lline2,NA,Lline2.0)
            Sline1 <- c(Sline1,NA,Sline1.0)
            Sline2 <- c(Sline2,NA,Sline2.0)
            Life_legend1 <- c(Life_legend1,rep(paste(c("RH at",unique(rbind(S,Suse)[,2])[i]," ",stressunit2),collapse = " "),101))
            Life_legend2 <- c(Life_legend2,rep(paste(c("Temperature at",unique(rbind(S,Suse)[,1])[i]," ",stressunit1),collapse = " "),101))
          }
        }
      }
      if(is.null(Suse)==TRUE){
        for(i in 1:length(unique(S[,1]))){
          if(dist=="Exponential"){
            Lline1.0<-lifestress.select(ls)[[1]](params,Sline1.0,unique(S[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params,unique(S[,1])[i],Sline2.0)
          }
          if(dist!="Exponential"){
            Lline1.0<-lifestress.select(ls)[[1]](params[2:length(params)],Sline1.0,unique(S[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params[2:length(params)],unique(S[,1])[i],Sline2.0)
          }
          if(i==1){
            Lline1 <- Lline1.0
            Lline2 <- Lline2.0
            Sline1 <- Sline1.0
            Sline2 <- Sline2.0
            Life_legend1 <- rep(paste(c("RH at",unique(S[,2])[i]," ",stressunit2),collapse = " "),100)
            Life_legend2 <- rep(paste(c("Temperature at",unique(S[,1])[i]," ",stressunit1),collapse = " "),100)
          } else{
            Lline1 <- c(Lline1,NA,Lline1.0)
            Lline2 <- c(Lline2,NA,Lline2.0)
            Sline1 <- c(Sline1,NA,Sline1.0)
            Sline2 <- c(Sline2,NA,Sline2.0)
            Life_legend1 <- c(Life_legend1,rep(paste(c("RH at",unique(S[,2])[i]," ",stressunit2),collapse = " "),101))
            Life_legend2 <- c(Life_legend2,rep(paste(c("Temperature at",unique(S[,1])[i]," ",stressunit1),collapse = " "),101))
          }
        }
      }
    }


    if(ls=="TempNonthermal" || ls=="Eyring3" || ls=="Eyring4"){
      # First plot line with changing temperature
      if(is.null(Suse)==FALSE){
        for(i in 1:length(unique(rbind(S,Suse)[,1]))){
          if(dist=="Exponential"){
            # unique(rbind(S,Suse)[,1])[i]
            # unique(rbind(S,Suse)[,2])[i]
            Lline1.0<-lifestress.select(ls)[[1]](params,Sline1.0,unique(rbind(S,Suse)[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params,unique(rbind(S,Suse)[,1])[i],Sline2.0)
          }
          if(dist!="Exponential"){
            Lline1.0<-lifestress.select(ls)[[1]](params[2:length(params)],Sline1.0,unique(rbind(S,Suse)[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params[2:length(params)],unique(rbind(S,Suse)[,1])[i],Sline2.0)
          }
          if(i==1){
            Lline1 <- Lline1.0
            Lline2 <- Lline2.0
            Sline1 <- Sline1.0
            Sline2 <- Sline2.0
            Life_legend1 <- rep(paste(c("Non-Thermal Stress at",unique(rbind(S,Suse)[,2])[i]," ",stressunit2),collapse = " "),100)
            Life_legend2 <- rep(paste(c("Temperature at",unique(rbind(S,Suse)[,1])[i]," ",stressunit1),collapse = " "),100)
          } else{
            Lline1 <- c(Lline1,NA,Lline1.0)
            Lline2 <- c(Lline2,NA,Lline2.0)
            Sline1 <- c(Sline1,NA,Sline1.0)
            Sline2 <- c(Sline2,NA,Sline2.0)
            Life_legend1 <- c(Life_legend1,rep(paste(c("Non-Thermal Stress at",unique(rbind(S,Suse)[,2])[i]," ",stressunit2),collapse = " "),101))
            Life_legend2 <- c(Life_legend2,rep(paste(c("Temperature at",unique(rbind(S,Suse)[,1])[i]," ",stressunit1),collapse = " "),101))
          }
        }
      }
      if(is.null(Suse)==TRUE){
        for(i in 1:length(unique(S[,1]))){
          if(dist=="Exponential"){
            Lline1.0<-lifestress.select(ls)[[1]](params,Sline1.0,unique(S[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params,unique(S[,1])[i],Sline2.0)
          }
          if(dist!="Exponential"){
            Lline1.0<-lifestress.select(ls)[[1]](params[2:length(params)],Sline1.0,unique(S[,2])[i])
            Lline2.0<-lifestress.select(ls)[[1]](params[2:length(params)],unique(S[,1])[i],Sline2.0)
          }
          if(i==1){
            Lline1 <- Lline1.0
            Lline2 <- Lline2.0
            Sline1 <- Sline1.0
            Sline2 <- Sline2.0
            Life_legend1 <- rep(paste(c("Non-Thermal Stress at",unique(S[,2])[i]," ",stressunit2),collapse = " "),100)
            Life_legend2 <- rep(paste(c("Temperature at",unique(S[,1])[i]," ",stressunit1),collapse = " "),100)
          } else{
            Lline1 <- c(Lline1,NA,Lline1.0)
            Lline2 <- c(Lline2,NA,Lline2.0)
            Sline1 <- c(Sline1,NA,Sline1.0)
            Sline2 <- c(Sline2,NA,Sline2.0)
            Life_legend1 <- c(Life_legend1,rep(paste(c("Non-Thermal Stress at",unique(S[,2])[i]," ",stressunit2),collapse = " "),101))
            Life_legend2 <- c(Life_legend2,rep(paste(c("Temperature at",unique(S[,1])[i]," ",stressunit1),collapse = " "),101))
          }
        }
      }
    }
    # Set up data frame for best fit line
    df_line <- data.frame(S = c(Sline1,Sline2), Sinv = c(1/Sline1,1/Sline2), L = c(Lline1,Lline2), best_fit = c(Life_legend1,Life_legend2), stress.group = c(rep(Slab,length(Sline1)),rep(Slab2,length(Sline2))))
  }
  # # Confidence upper and lower bound
  # # Mean square error
  # MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
  # # SS
  # SS <- sum((1)^2)
  # CONFDIFF <- qt(confid,dim(data)[1])*sqrt(MSE*(1/dim(data)[1] + (Lline - mean(Lline))^2))
  # df_confid_bound <- data.frame(S = Sline, Sinv = 1/Sline, Llower = Lline, Lupper = Lline, best_fit = rep("Fitted",100))
  # # Predicative upper and lower bound
  # df <- data.frame(X = xrange, YCDF = ycdf, YCDFlow = ycdf_low, YCDFhigh = ycdf_high, best_fit = rep("Fitted",1000))

  # UPDATE (11/9/2023): Adding the test of the plots now.  These will be standard output.  Going to test Arrhenius first
  # for all conditions and then do the rest when I'm satisfied.

  if (ls=="Linear") {
    # theta[1] - parameter a, theta[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = Lline - CONFDIFF, Lupper = Lline + CONFDIFF, best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = Lline - PREDICTDIFF, Lupper = Lline + PREDICTDIFF, best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="Exponential"){
    # theta[1] - parameter a, theta[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="Exponential2"){
    # theta[1] - parameter a, theta[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - 1/S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    # SS
    SS <- sum(((1/sort.xircstressdata(data)[[3]]) - mean(1/sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + ((((1/Sline) - mean(1/(sort.xircstressdata(data)[[3]]))))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Sinv = (1/Sline), Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + ((((1/Sline) - mean(1/(sort.xircstressdata(data)[[3]]))))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Sinv = (1/Sline), Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.Sinv, xmax = c(x.shift.Sinv+dist.plot.y*0.8*(dist_scale_invS/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }

  }

  if (ls=="Arrhenius") {
    # lsparams[1] - parameter Ea, lsparams[2] - parameter b
    # Temperature HaS to be in Kelvin for this to work
    K<-8.617385e-5
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - 1/S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    # SS
    SS <- sum(((1/sort.xircstressdata(data)[[3]]) - mean(1/sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + ((((1/Sline) - mean(1/(sort.xircstressdata(data)[[3]]))))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Sinv = (1/Sline), Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + ((((1/Sline) - mean(1/(sort.xircstressdata(data)[[3]]))))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Sinv = (1/Sline), Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.Sinv, xmax = c(x.shift.Sinv+dist.plot.y*0.8*(dist_scale_invS/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="Eyring") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - 1/S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.Sinv, xmax = c(x.shift.Sinv+dist.plot.y*0.8*(dist_scale_invS/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="Eyring2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - 1/S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.Sinv, xmax = c(x.shift.Sinv+dist.plot.y*0.8*(dist_scale_invS/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=Sinv,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=Sinv, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=Sinv, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="Power") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="PowerwithBias") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = Lline - CONFDIFF, Lupper = Lline + CONFDIFF, best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = Lline - PREDICTDIFF, Lupper = Lline + PREDICTDIFF, best_fit = rep("Confidence",100))
    }

    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="InversePower") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="InversePower2") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    # return(list(x.shift.S,dist.plot.y,dist_legend))
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0)) +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }
  if (ls=="Logarithmic") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    if(is.null(Slab)==TRUE){
      Slab <- "Characteristic Stress - S"
    }
    if(is.null(Llab)==TRUE){
      Llab <- "Characteristic Life - L"
    }

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    if(dist=="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
    }
    if(dist!="Exponential"){
      MSE <- (1/(dim(data)[1]-2))*sum((lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]]) - sort.xircstressdata(data)[[1]])^2)
    }
    # SS
    SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    df_confid_bound <- data.frame(S = Sline, Llower = Lline - CONFDIFF, Lupper = Lline + CONFDIFF, best_fit = rep("Confidence",100))
    # Predicative upper and lower bound
    if(is.null(predic_int)==FALSE){
      PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
      df_predict_bound <- data.frame(S = Sline, Llower = Lline - PREDICTDIFF, Lupper = Lline + PREDICTDIFF, best_fit = rep("Confidence",100))
    }
    # PDF distribution overlay by stress (NEW)
    df_PDF_imposed <- data.frame(xmin = x.shift.S, xmax = c(x.shift.S+dist.plot.y*0.8*(dist_scale_S/max(dist.plot.y))), y = dist.plot.x, PDFfit = dist_legend)

    relationplot<-ggplot() +
      geom_line(data=df_line, aes(x=S,y=L), colour = 'black', linewidth = 0.9, linetype = "dashed") +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      scale_x_continuous(expand=c(0, 0),trans = 'log10') +
      scale_y_continuous(expand=c(0, 0)) +
      xlab(Slab) +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }
  if (ls=="TempHumidity") {
    # lsparams[1] - parameter A, lsparams[2] - parameter a, lsparams[3] - parameter b

    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    # if(dist=="Exponential"){
    #   MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    # }
    # if(dist!="Exponential"){
    #   MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    # }
    # # SS
    # SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    # CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    # df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # # Predicative upper and lower bound
    # if(is.null(predic_int)==FALSE){
    #   PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    #   df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    # }
    # PDF distribution overlay by stress (NEW)

    df_PDF_imposed <- data.frame(xmin = c(x.shift.S1,x.shift.S2),
                                 xmax = c(x.shift.S1 + dist.plot.y*0.8*(dist_scale_S1/max(dist.plot.y)),
                                          x.shift.S2 + dist.plot.y*0.8*(dist_scale_S2/max(dist.plot.y))),
                                 y = c(dist.plot.x,dist.plot.x), PDFfit = c(dist_legend,dist_legend),
                                 stress.group = c(rep(Slab,length(dist_legend)),rep(Slab2,length(dist_legend))))

    # posterior_beta <- density(c(fit$draws("beta")))
    # posterior_alpha <- density(c(fit$draws("alpha")))
    # df_posterior <- data.frame(x = c(posterior_alpha$x,posterior_beta$x),
    #                            ymin = rep(0,(length(posterior_alpha$x)+length(posterior_beta$x))),
    #                            ymax = c(posterior_alpha$y,posterior_beta$y),
    #                            distlabel = c(rep("α",length(posterior_alpha$x)),rep("β",length(posterior_beta$x))))
    #
    # # Density plot for Exponential rate parameter posterior
    # plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
    #   theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
    #   facet_wrap(~distlabel, dir="v", scales = "free") +
    #   scale_x_continuous(expand=c(0, 0)) +
    #   scale_y_continuous(expand=c(0, 0)) +
    #   xlab(" ")
    # ylab("density")

    relationplot<-ggplot() +
      geom_path(data=df_line, aes(x=S,y=L, linetype = best_fit), colour = 'black', linewidth = 0.9) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~stress.group, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0),trans = 'log10') +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(" ") +
      ylab(Llab)

    # return(relationplot)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        facet_wrap(~stress.group, dir="v", scales = "free") +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        facet_wrap(~stress.group, dir="v", scales = "free") +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if (ls=="TempNonthermal") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b, lsparams[3] - parameter c


    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    # if(dist=="Exponential"){
    #   MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    # }
    # if(dist!="Exponential"){
    #   MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    # }
    # # SS
    # SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    # CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    # df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # # Predicative upper and lower bound
    # if(is.null(predic_int)==FALSE){
    #   PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    #   df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    # }
    # PDF distribution overlay by stress (NEW)

    df_PDF_imposed <- data.frame(xmin = c(x.shift.S1,x.shift.S2),
                                 xmax = c(x.shift.S1 + dist.plot.y*0.8*(dist_scale_S1/max(dist.plot.y)),
                                          x.shift.S2 + dist.plot.y*0.8*(dist_scale_S2/max(dist.plot.y))),
                                 y = c(dist.plot.x,dist.plot.x), PDFfit = c(dist_legend,dist_legend),
                                 stress.group = c(rep(Slab,length(dist_legend)),rep(Slab2,length(dist_legend))))

    relationplot<-ggplot() +
      geom_path(data=df_line, aes(x=S,y=L, linetype = best_fit), colour = 'black', linewidth = 0.9) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~stress.group, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0),trans = 'log10') +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(" ") +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        facet_wrap(~stress.group, dir="v", scales = "free") +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        facet_wrap(~stress.group, dir="v", scales = "free") +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }
  if (ls=="Eyring3") {
    # lsparams[1] - parameter a, lsparams[2] - parameter b
    # lsparams[3] - parameter c, lsparams[4] - parameter d


    # Confidence upper and lower bound
    # Mean square error
    # y = ln(L)
    # if(dist=="Exponential"){
    #   MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params,sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    # }
    # if(dist!="Exponential"){
    #   MSE <- (1/(dim(data)[1]-2))*sum((log(lifestress.select(ls)[[1]](params[2:length(params)],sort.xircstressdata(data)[[3]])) - log(sort.xircstressdata(data)[[1]]))^2)
    # }
    # # SS
    # SS <- sum(((sort.xircstressdata(data)[[3]]) - mean(sort.xircstressdata(data)[[3]]))^2)
    # CONFDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*((1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    # df_confid_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - CONFDIFF), Lupper = exp(log(Lline) + CONFDIFF), best_fit = rep("Confidence",100))
    # # Predicative upper and lower bound
    # if(is.null(predic_int)==FALSE){
    #   PREDICTDIFF <- qt(confid,(dim(data)[1]-2))*sqrt(MSE*(1 + (1/dim(data)[1]) + (((Sline - mean(sort.xircstressdata(data)[[3]])))^2)/SS))
    #   df_predict_bound <- data.frame(S = Sline, Llower = exp(log(Lline) - PREDICTDIFF), Lupper = exp(log(Lline) + PREDICTDIFF), best_fit = rep("Confidence",100))
    # }
    # PDF distribution overlay by stress (NEW)

    df_PDF_imposed <- data.frame(xmin = c(x.shift.S1,x.shift.S2),
                                 xmax = c(x.shift.S1 + dist.plot.y*0.8*(dist_scale_S1/max(dist.plot.y)),
                                          x.shift.S2 + dist.plot.y*0.8*(dist_scale_S2/max(dist.plot.y))),
                                 y = c(dist.plot.x,dist.plot.x), PDFfit = c(dist_legend,dist_legend),
                                 stress.group = c(rep(Slab,length(dist_legend)),rep(Slab2,length(dist_legend))))

    relationplot<-ggplot() +
      geom_path(data=df_line, aes(x=S,y=L, linetype = best_fit), colour = 'black', linewidth = 0.9) +
      theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
      facet_wrap(~stress.group, dir="v", scales = "free") +
      scale_x_continuous(expand=c(0, 0),trans = 'log10') +
      scale_y_continuous(expand=c(0, 0),trans = 'log10') +
      xlab(" ") +
      ylab(Llab)

    if(is.null(confid_int)==FALSE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }
    if(is.null(confid_int)==TRUE && is.null(predic_int)==FALSE){
      relationplot <- relationplot + geom_ribbon(data = df_predict_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "red",alpha = 0.25)
    }
    if(is.null(confid_int)==FALSE && is.null(predic_int)==TRUE){
      relationplot <- relationplot + geom_ribbon(data = df_confid_bound, aes(x=S, ymin = Llower, ymax = Lupper), fill = "blue",alpha = 0.25)
    }

    if(is.null(Suse) == FALSE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)],11), breaks = data_breaks) +
        scale_color_manual(values=c(rep(col_legend[1:length(fulldatabystress)],2),"black"))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:(1+length(fulldatabystress))], breaks = pdf_breaks)
    }
    if(is.null(Suse) == TRUE){
      # Add data points
      relationplot <- relationplot + geom_point(data=df_data, aes(x=S, y=L, shape=data.group), color="black", size=3) +
        facet_wrap(~stress.group, dir="v", scales = "free") +
        scale_shape_manual("Data Points",values=c(shape_legend[1:length(fulldatabystress)],shape_legend2[1:length(fulldatabystress)]), breaks = data_breaks) +
        scale_color_manual(values=rep(col_legend[1:length(fulldatabystress)],2))
      # Add PDF distribution bands
      relationplot <- relationplot + geom_ribbon(data=df_PDF_imposed, aes(xmin = xmin,xmax = xmax, y = y, fill = PDFfit), alpha=0.25) +
        facet_wrap(~stress.group, dir="v", scales = "free") +
        scale_fill_manual(paste(c(dist,"PDF"),collapse = " "), values=col_legend[1:length(fulldatabystress)], breaks = pdf_breaks)
    }
  }

  if(is.null(Suse) == FALSE){
    # return(list(df_data,df_line,relationplot=relationplot,Luse=Luse))
    return(list(relationplot=relationplot,Luse=Luse))
  }
  if(is.null(Suse) == TRUE){
    # return(list(df_data,df_line,relationplot=relationplot))
    return(list(relationplot=relationplot))
  }

}
