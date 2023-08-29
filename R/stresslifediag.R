# Stress-Life Diagram
# Developed by Reuel Smith, 2022

SN.diagram <- function(input_type,dat,stressunits,options){
  # Define input type for data: 1 - data points, 2- cyclic stress range, 3 - multi-axial stress
  # Define units in either MPa (default) or psi.  Database will have units for both.
  # Options for data type 1 to include: Su, Se, BHN, fatigue strength
  library(pracma)
  library(ggplot2)
  library(plyr)
  suppressWarnings({
  # Check units and set up axis labels
  if(missing(stressunits) || stressunits == 1){
    stressunits <- c("MPa")
    Su_threshold <- 200*6.8947572932
    Se_threshold <- 100*6.8947572932
    BHN_mult <- 0.25*6.8947572932
  }
  if(stressunits == 2){
    stressunits <- c("ksi")
    Su_threshold <- 200
    Se_threshold <- 100
    BHN_mult <- 0.25
  }

  # Initialize Su, Se, and BHN
  Se <- logical(0)
  Su <- logical(0)
  BHN <- logical(0)

  # Initialize A and b
  A <- logical(0)
  b <- logical(0)

  # Initialize Stress and Life
  S_ <- logical(0)
  N_ <- logical(0)
  runoff_S <- logical(0)
  runoff_N <- logical(0)

  # Initialize Stress Trace and Life Trace
  S_trace <- logical(0)
  Sar <- logical(0)
  N_trace <- logical(0)
  Scurve_trace1 <- logical(0)
  Ncurve_trace1 <- logical(0)
  Scurve_trace2 <- logical(0)
  Ncurve_trace2 <- logical(0)

  # Check units and set up axis labels
  Xlab <- paste(c("Life to Failure, N (Cycles)"),collapse="", sep = "_")
  # Ylab <- bquote(.(c("Alternating Stress, S[a]", "(",stressunits,")")))
  Ylab <- paste(c("Alternating Stress, S (",stressunits,")"),collapse="", sep = "_")

  # ==================================
  # Calculation Type 1 - Data Points
  # ==================================
  if(input_type==1){
    # Check to see if data is in list form (if it isn't NA that is)
    if(isFALSE(is.na(dat)[1]) && is.list(dat)==FALSE){
      stop('Enter S-N data as a list of vectors.')
    }
  }

  # ==================================
  # Calculation Type 2 - Cyclic Stress Range OR
  # Calculation Type 3 - Multi-Axial Stress
  # ==================================
  # Both can be used for fully reversed cycles, but when not the case,
  # need a model to base estimate upon
  if(input_type==2 || input_type==3){
    # Compute the stress amplitude and the mean stress based on input_type
    # ==================================
    # TYPE 2
    # ==================================
    if(input_type==2){
      if(missing(options)==FALSE){
        # Pull Smax and Smin if available
        if(length(options$Srange)>0){
          if(is.na(options$Srange)[1]==0){
            Smax <- max(options$Srange)
            Smin <- min(options$Srange)
          }
        } else {
          stop('Enter options for Srange (c(Smin,Smax)).')
        }
      } else {
        stop('Enter options for Srange (c(Smin,Smax)).')
      }
      # Stress Amplitude
      Sa <- (Smax - Smin)*0.5
      # Mean Stress
      Sm <- (Smax + Smin)*0.5
    }
    # ==================================
    # TYPE 3
    # ==================================
    if(input_type==3){
      # Initialize multi-axial stresses
      sig_xa <- 0
      sig_ya <- 0
      sig_za <- 0
      tau_xya <- 0
      tau_yza <- 0
      tau_xza <- 0
      sig_xm <- 0
      sig_ym <- 0
      sig_zm <- 0
      # ==================================
      # Find and replace any given values from options
      # ==================================
      if(isFALSE(missing(options))){
        # NORMAL STRESS AMPLITUDE
        if(length(options$sig_xa)==1){
          sig_xa <- options$sig_xa
        }
        if(length(options$sig_ya)==1){
          sig_ya <- options$sig_ya
        }
        if(length(options$sig_za)==1){
          sig_za <- options$sig_za
        }
        # NORMAL MEAN STRESS
        if(length(options$sig_xm)==1){
          sig_xm <- options$sig_xm
        }
        if(length(options$sig_ym)==1){
          sig_ym <- options$sig_ym
        }
        if(length(options$sig_zm)==1){
          sig_zm <- options$sig_zm
        }
        # SHEAR STRESS AMPLITUDE
        if(length(options$tau_xya)==1){
          tau_xya <- options$tau_xya
        }
        if(length(options$tau_yza)==1){
          tau_yza <- options$tau_yza
        }
        if(length(options$tau_xza)==1){
          tau_xza <- options$tau_xza
        }
        # ==================================
        # You have an option however to enter the axial stress ranges if known
        # Then the amplitude and mean stresses will be computed as normal
        # ==================================
        # NORMAL STRESS AMPLITUDE AND MEAN
        if(length(options$Srangex) == 2){
          sig_xa <- (max(options$Srangex) - min(options$Srangex))*0.5
          sig_xm <- (max(options$Srangex) + min(options$Srangex))*0.5
        }
        if(length(options$Srangey) == 2){
          sig_ya <- (max(options$Srangey) - min(options$Srangey))*0.5
          sig_ym <- (max(options$Srangey) + min(options$Srangey))*0.5
        }
        if(length(options$Srangez) == 2){
          sig_za <- (max(options$Srangez) - min(options$Srangez))*0.5
          sig_zm <- (max(options$Srangez) + min(options$Srangez))*0.5
        }
        # SHEAR STRESS AMPLITUDE
        if(length(options$Trangexy) == 2){
          tau_xya <- (max(options$Trangexy) - min(options$Trangexy))*0.5
        }
        if(length(options$Trangeyz) == 2){
          tau_yza <- (max(options$Trangeyz) - min(options$Trangeyz))*0.5
        }
        if(length(options$Trangexz) == 2){
          tau_xza <- (max(options$Trangexz) - min(options$Trangexz))*0.5
        }
        # Check to see if anything was entered at all
        if(sum(sig_xa,sig_ya,sig_za,tau_xya,tau_yza,tau_xza,sig_xm,sig_ym,sig_zm) == 0){
          stop('Enter options for axial stress amplitudes, means, and shear stress amplitude or their loading ranges.')
        }
      }
      # Stress Amplitude
      Sa <- (2^-0.5)*sqrt(sum((sig_xa-sig_ya)^2,(sig_ya-sig_za)^2,(sig_xa-sig_za)^2,6*((sum(tau_xya,tau_yza,tau_xza))^2)))
      # Mean Stress
      Sm <- sum(sig_xm,sig_ym,sig_zm)
    }
  }

  if(is.na(dat)[1]==FALSE){
    # Set input data to SN data
    SNdat <- dat

    # Pull S and N data if available
    S_ <- SNdat[[1]]
    N_ <- SNdat[[2]]
    # Check for and pull runoff data
    if(length(dat)==4){
      runoff_S <- dat[[3]]
      runoff_N <- dat[[4]]
    }
  }

  # Pull A, b, Se, Su, and/or BHN if available
  if(isFALSE(missing(options))){
    if(length(options$A)==1){
      A <- options$A
    }
    if(length(options$b)==1){
      b <- options$b
    }
    # NOTE: A and b will only be used if no data is given as input
    if(length(options$Se)==1){
      Se <- options$Se
    }
    if(length(options$Su)==1){
      Su <- options$Su
    }
    if(length(options$BHN)==1){
      BHN <- options$BHN
    }
  }

  # STREAMLINE PART 2 (9/6/2022)
  # Find Endurance Limit Se
  if(length(Se)<=1 && length(runoff_S)>=2){
    # By mean or average of runoff stress data if it exists
    # First pass is meant to override the entered Se if it exists
    Se <- mean(runoff_S)
  }
  if(length(Se)==1 && length(runoff_S)<2){
    # SET to Se already so no computation needed
  }
  if(length(Se)==0 && length(Su)==1){
    # Compute estimate for endurance limit based on ultimate stress Su
    if(Su <= Su_threshold){
      Se <- 0.5*Su
    } else{
      Se <- Se_threshold
    }
  }
  if(length(Se)==0 && length(BHN)==1){
    # Compute estimate for endurance limit based on BHN
    if(BHN <= 400){
      Se <- BHN_mult*BHN
    } else{
      Se <- Se_threshold
    }
  }

  # Compute A and b
  if(is.na(dat)[1]==FALSE){
    if(length(S_)==1){
      # If data only has one point, assume 10^6 as your N_Se
      S_ <- c(S_,Se)
      N_ <- c(N_,10^6)
    }
    if(length(S_)>=2){
      # Pull A and b from data if available
    }
    logS <- log(S_)
    logN <- log(N_)
    # Compute A and b for S_a = A N^b
    params <- lm(logS ~ poly(logN, 1, raw=TRUE))
    A <- exp(summary(params)$coefficients[1,1])
    b <- summary(params)$coefficients[2,1]
  }


  if(length(S_)==0 && sum(length(A),length(b))==2){
    # A and b given will be used
    # Final A and b condition where one data is given to follow
  }
  # # Check and see if we have an A and b
  # if(length(S_)==0 && sum(length(A),length(b)) < 2){
  #   stop('Enter either stress and life data in dat or enter an A and b variable set in options.')
  # }

  # Then if no runoff data exists, retroactively insert some
  if(length(runoff_S)==0){
    runoff_S <- Se
    runoff_N <- 10^7
  }

  # Compute S_1000 only if Su is not known
  S1000 <- A*(1000^b)
  if(length(Su)==0){
    Su <-S1000/0.9
    # Retroactively compute Se if it still does not exist
    if(length(Se)==0){
      if(Su <= Su_threshold){
        Se <- 0.5*Su
      } else{
        Se <- Se_threshold
      }
    }
  }

  # ==================================
  # Second call for input_type 2 and 3 to obtain the equivalent stress sig_ar
  if(input_type==2 || input_type==3){
    if(isFALSE(missing(options))){
      # Check correction relationships options$corr_rel and needed input
      if(length(is.na(options$corr_rel))==0){
        stop('Please enter a correction relationship as "corr_rel".')
      } else{
        if(options$corr_rel == "Soderberg" && length(options$Sy) == 1){
          Sy <- options$Sy
          Sar <- Sa/(1 - (Sm/Sy))
        }
        if(options$corr_rel == "ModifiedGoodman" && length(Su) == 1){
          Sar <- Sa/(1 - (Sm/Su))
        }
        if(options$corr_rel == "Morrow" && length(options$sigf) == 1){
          sigf <- options$sigf
          Sar <- Sa/(1 - (Sm/sigf))
        }
        if(options$corr_rel == "Gerber" && length(Su) == 1){
          Sar <- Sa/(1 - ((Sm/Su)^2))
        }
        if(options$corr_rel == "SWT"){
          # Stress Ratio
          R <- Smin/Smax
          Sar <- Sa/sqrt(0.5*(1 - R))
        }
        if(options$corr_rel == "Walker" && length(options$gam) == 1){
          # Stress Ratio
          R <- Smin/Smax
          Sar <- Sa/((0.5*(1 - R))^(1-gam))
        }
      }
      }
  }

  N_Se <- (Se/A)^(1/b)
  Ncurve <- c(linspace(1000,N_Se,100),max(runoff_N))
  Scurve <-c(A*((linspace(1000,N_Se,100))^b),Se)

  if(input_type==1){
    # Establish the stress and life trace if given
    if(isFALSE(missing(options))){
      if(length(options$Strace)==1){
        S_trace <- options$Strace
        N_Sar1 <- (S_trace/A)^(1/b)
        Ncurve_trace1 <- c(1000,N_Sar1,N_Sar1)
        Scurve_trace1 <- c(S_trace,S_trace,round_any(Se, 10^floor(log10(Se)), f = floor))
      }
      if(length(options$Ntrace)==1){
        N_trace <- options$Ntrace
        Sar1 <- A*(N_trace^b)
        Ncurve_trace2 <- c(1000,N_trace,N_trace)
        Scurve_trace2 <- c(Sar1,Sar1,round_any(Se, 10^floor(log10(Se)), f = floor))
      }
    }
  }
  if(input_type==2 || input_type==3){
    # Establish the stress and life trace if given
    if(isFALSE(missing(options))){
      N_Sar1 <- (Sar/A)^(1/b)
      S_trace <- Sar
      Ncurve_trace1 <- c(1000,N_Sar1,N_Sar1)
      Scurve_trace1 <- c(Sar,Sar,round_any(Se, 10^floor(log10(Se)), f = floor))
      if(length(options$Ntrace)==1){
        N_trace <- options$Ntrace
        Sar1 <- A*(N_trace^b)
        Ncurve_trace2 <- c(1000,N_trace,N_trace)
        Scurve_trace2 <- c(Sar1,Sar1,round_any(Se, 10^floor(log10(Se)), f = floor))
      }
    }
  }

  if(length(S_)>2){
    # Compute upper and lower bounds
    res_logS <- logS - (log(A) + b*logN)
    logCI_S <- quantile(res_logS, 0.995)
    Scurvelow <- exp(c(log(A) + b*log((linspace(1000,N_Se,100))) - logCI_S,log(A) + b*log(N_Se) - logCI_S))
    Scurvehigh <- exp(c(log(A) + b*log((linspace(1000,N_Se,100))) + logCI_S,log(A) + b*log(N_Se) + logCI_S))

    if(sum(length(Ncurve_trace1),length(Ncurve_trace2))==0){
      df<-data.frame(S=c(S_,rep(NA,length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve), Slineupper = c(rep(NA,length(S_)+length(runoff_N)),Scurvehigh), Slinelower = c(rep(NA,length(S_)+length(runoff_N)),Scurvelow), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve))))
    }
    if(length(Ncurve_trace1)>0 && length(Ncurve_trace2)>0){
      df<-data.frame(S=c(S_,rep(NA,6+length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,6+length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,6+length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,6+length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve,rep(NA,6)), Slineupper = c(rep(NA,length(S_)+length(runoff_N)),Scurvehigh,rep(NA,6)), Slinelower = c(rep(NA,length(S_)+length(runoff_N)),Scurvelow,rep(NA,6)), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve,rep(NA,6)), Straceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace1,rep(NA,3)), Ntraceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace1,rep(NA,3)), Straceline2 = c(rep(NA,3+length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace2), Ntraceline2 = c(rep(NA,3+length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace2), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,6+length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve)),rep(NA,6)))
    }
    if(length(Ncurve_trace1)>0 && length(Ncurve_trace2)==0){
      df<-data.frame(S=c(S_,rep(NA,3+length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,3+length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,3+length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,3+length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve,rep(NA,3)), Slineupper = c(rep(NA,length(S_)+length(runoff_N)),Scurvehigh,rep(NA,3)), Slinelower = c(rep(NA,length(S_)+length(runoff_N)),Scurvelow,rep(NA,3)), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve,rep(NA,3)), Straceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace1), Ntraceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace1), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,3+length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve)),rep(NA,3)))
    }
    if(length(Ncurve_trace1)==0 && length(Ncurve_trace2)>0){
      df<-data.frame(S=c(S_,rep(NA,3+length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,3+length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,3+length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,3+length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve,rep(NA,3)), Slineupper = c(rep(NA,length(S_)+length(runoff_N)),Scurvehigh,rep(NA,3)), Slinelower = c(rep(NA,length(S_)+length(runoff_N)),Scurvelow,rep(NA,3)), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve,rep(NA,3)), Straceline2 = c(rep(NA,length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace2), Ntraceline2 = c(rep(NA,length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace2), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,3+length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve)),rep(NA,3)))
    }
  } else {
    # No upper and lower bounds for two points
    if(sum(length(Ncurve_trace1),length(Ncurve_trace2))==0){
      df<-data.frame(S=c(S_,rep(NA,length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve))))
    }
    if(length(Ncurve_trace1)>0 && length(Ncurve_trace2)>0){
      df<-data.frame(S=c(S_,rep(NA,6+length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,6+length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,6+length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,6+length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve,rep(NA,6)), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve,rep(NA,6)), Straceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace1,rep(NA,3)), Ntraceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace1,rep(NA,3)), Straceline2 = c(rep(NA,3+length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace2), Ntraceline2 = c(rep(NA,3+length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace2), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,6+length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve)),rep(NA,6)))
    }
    if(length(Ncurve_trace1)>0 && length(Ncurve_trace2)==0){
      df<-data.frame(S=c(S_,rep(NA,3+length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,3+length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,3+length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,3+length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve,rep(NA,3)), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve,rep(NA,3)), Straceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace1), Ntraceline1 = c(rep(NA,length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace1), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,3+length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve)),rep(NA,3)))
    }
    if(length(Ncurve_trace1)==0 && length(Ncurve_trace2)>0){
      df<-data.frame(S=c(S_,rep(NA,3+length(runoff_S)+length(Ncurve))), N = c(N_,rep(NA,3+length(runoff_S)+length(Ncurve))), Srunoff=c(rep(NA,length(S_)),runoff_S,rep(NA,3+length(Ncurve))), Nrunoff = c(rep(NA,length(N_)),runoff_N,rep(NA,3+length(Ncurve))), Sline = c(rep(NA,length(S_)+length(runoff_N)),Scurve,rep(NA,3)), Nline = c(rep(NA,length(S_)+length(runoff_N)),Ncurve,rep(NA,3)), Straceline2 = c(rep(NA,length(S_)+length(runoff_N)+length(Scurve)),Scurve_trace2), Ntraceline2 = c(rep(NA,length(S_)+length(runoff_N)+length(Ncurve)),Ncurve_trace2), data_points = c(rep("failed",length(S_)),rep("survived",length(runoff_N)),rep(NA,3+length(Ncurve))), datapt_or_curvefit = c(rep(NA,length(S_)+length(runoff_N)),rep("curvefit",length(Ncurve)),rep(NA,3)))
    }
  }

  plotout<-ggplot() +
    geom_point(data=df, aes(N,S), colour = 'red', size = 1.9) +
    geom_point(data=df, aes(Nrunoff,Srunoff), colour = 'green4', shape=17, size = 1.9) +
    geom_line(data=df, aes(Nline,Sline), colour = "black", size = 0.9) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab(Xlab) +
    ylab(Ylab)

  # Plot boundary line if data is greater than 2
  if(length(S_)>2){
    plotout<-plotout + geom_ribbon(data=df, aes(ymin = Slinelower,ymax = Slineupper,x=Nline), alpha=0.25,fill = "red")
  }
  # Plot stress trace if provided
  if(length(Ncurve_trace1)>0){
    plotout<-plotout + geom_line(data=df, aes(Ntraceline1,Straceline1), colour = "red", size = 0.9, linetype = "dashed")
  }
  # Plot cycle trace if provided
  if(length(Ncurve_trace2)>0){
    plotout<-plotout + geom_line(data=df, aes(Ntraceline2,Straceline2), colour = "blue", size = 0.9, linetype = "dashed")
  }

  # Produce some output text that summarizes the results
  if(is.na(dat)[1]==FALSE){
    cat(c("The estimate for S-N curve S = ",A," x N_f^",b,".\n"),sep = "")
  }
  if(is.na(dat)[1]==TRUE && length(options$A)==1 && length(options$b)==1){
    cat(c("The given S-N curve is S = ",A," x N_f^",b,".\n"),sep = "")
  }
  if(length(runoff_S)>=2 || (length(options$Su)==1 && length(options$Se)<1)){
    cat(c("The estimate for endurance limit is S_e = ",Se," ",stressunits,".\n\n"),sep = "")
  }
  if(isFALSE(missing(options))){
    if(length(options$Se)==1){
      cat(c("The given endurance limit is S_e = ",Se," ",stressunits,".\n\n"),sep = "")
    }

    if((length(options$Strace)==1 && options$Strace > Se) || (input_type==2 || input_type==3)){
      cat(c("The life estimate at ",S_trace," ",stressunits," is ",N_Sar1," cycles.\n"),sep = "")
    }
    if((length(options$Strace)==1 && options$Strace <= Se)){
      cat(c("Enter a trace stress greater than the endurance limit S_e = ",Se," ",stressunits,".\n"),sep = "")
    }
    if(length(options$Ntrace)==1){
      cat(c("The alternating stress estimate at ",N_trace," cycles is ",Sar1," ",stressunits,".\n\n"),sep = "")
    }
  }


  if(sum(length(S_trace),length(N_trace),length(Sar))==0){
    return(list(SNdiag = plotout,A = A,b = b,Se = Se,Su = Su))
  }

  if(length(options$Ntrace)==1){
    return(list(SNdiag = plotout,A = A,b = b,Se = Se,Su = Su, Ntrace = N_trace, Sequiv = Sar1))
  }
  if((length(options$Strace)==1 && options$Strace > Se) || (input_type==2 || input_type==3)){
    return(list(SNdiag = plotout,A = A,b = b,Se = Se,Su = Su, Strace = S_trace, Nequiv = N_Sar1))
  }
  if(length(options$Ntrace)==1 && ((length(options$Strace)==1 && options$Strace > Se) || (input_type==2 || input_type==3))){
    return(list(SNdiag = plotout,A = A,b = b,Se = Se,Su = Su, Ntrace = N_trace, Sequiv = Sar1, Strace = S_trace, Nequiv = N_Sar1))
  }
  })
}
