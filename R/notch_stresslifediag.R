# Stress-Life Diagram with Notch Effect
# Developed by Reuel Smith, 2022

notch.SN.diagram <- function(dimensions,geometry,stressunits,loadconditions){
  # from example 2.7
  # dimensions = list(W = 100 (mm), r = 10 (mm), t = 5 (mm)), geometry = "rect_2semicirc_edge",
  # P or Smax = 150000 (N), Su = 790 (MPa), a = 0.145 (mm), ratioKf = 0.3
  library(pracma)
  library(ggplot2)
  library(plyr)

  if(length(loadconditions$Su)==0){
    stop('Please enter an ultimate strength, Su.')
  }

  # Check units and set up axis labels
  if(missing(stressunits) || stressunits == 1){
    Su_threshold <- 200*6.8947572932
    Se_threshold <- 100*6.8947572932
    BHN_mult <- 0.25*6.8947572932
    stressunitslabel <- "MPa"
  }
  if(stressunits == 2){
    Su_threshold <- 200
    Se_threshold <- 100
    BHN_mult <- 0.25
    stressunitslabel <- "ksi"
  }

  # Check units and set up axis labels
  Xlab <- paste(c("Life to Failure, N (Cycles)"),collapse="", sep = "_")
  Ylab <- paste(c("Net Alternating Stress, S_a (",stressunitslabel,")"),collapse="", sep = "_")

  # Compute the stress concentration factor, Kt
  Kt <- stress.concentration.factor(dimensions,geometry)

  # Compute the notch sensitivity factor, q
  # metric (mm)
  if(length(loadconditions$a)!=1 && stressunits == 1){
    a <- 0.0254*(((300*6.8947572932)/loadconditions$Su)^1.8)
  }
  # English (inches)
  if(length(loadconditions$a)==0 && stressunits == 2){
    a <- (1e-3)*((300/loadconditions$Su)^1.8)
  }
  if(length(loadconditions$a)==1){
    a <- loadconditions$a
  }
  q <- 1/(1 + (a/dimensions$r))

  # Compute fatigue stress concentration factor, Kf
  Kf <- 1 + q*(Kt - 1)

  # Compute Kf ratio
  Kfratio <-fatigue.stress.concentration.factor.ratio(loadconditions$Su,"Steel",stressunits)

  # Compute fatigue stress concentration factor at 1000 cycles, Kf'
  # Kfp <- (Kf - 1)*loadconditions$Kfratio + 1
  Kfp <- (Kf - 1)*Kfratio + 1

  # Find endurance limit unless it is given
  if(length(loadconditions$Se) == 1){
    Se <- loadconditions$Se
  }
  if(length(loadconditions$Se) == 0 && length(loadconditions$Su) == 1 && loadconditions$Su <= Su_threshold){
    Se <- 0.5*loadconditions$Su
  }
  if(length(loadconditions$Se) == 0 && length(loadconditions$Su) == 1 && loadconditions$Su > Su_threshold){
    Se <- Se_threshold
  }

  # Approximate S_1000 as a function of Su
  S1000 <- 0.9*loadconditions$Su

  # Compute notched Se, S1, and S1000
  Se_notch <- Se/Kf
  S1000_notch <- S1000/Kfp

  # Metric
  if(stressunits == 1){
    S1_notch <- loadconditions$Su + 50*6.8947572932
  }
  # English
  if(stressunits == 2){
    S1_notch <- loadconditions$Su + 50
  }
  # NOTE: 50 ksi is the approximation corrector to attain true fracture strength for STEELS only.
  # Future iterations of this code will need to make allowance for material input so that this may
  # be adjusted accordingly

  # Calculate Life to Failure based on sample type and area
  if(geometry == "rect_1semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$t) == 1 && length(dimensions$r) == 1){
    Area <- dimensions$t*(dimensions$W - dimensions$r)
  }
  if(geometry == "rect_2semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$t) == 1 && length(dimensions$r) == 1){
    Area <- dimensions$t*(dimensions$W - 2*dimensions$r)
  }
  Snet <- loadconditions$Sa/Area

  # ==================================
  # Plotting
  # ==================================
  Ncurve <- c(1,1000,10^6,10^9)
  Scurve_unnotched <- c(S1_notch,S1000,Se,Se)
  Scurve_notched <- c(S1_notch,S1000_notch,Se_notch,Se_notch)
  if(Snet > S1000_notch && Snet <= S1_notch){
    logN <- log(c(1,1000))
    logS <- log(c(S1_notch,S1000_notch))
  }
  if(Snet >= Se_notch && Snet <= S1000_notch){
    logN <- log(c(1000,10^6))
    logS <- log(c(S1000_notch,Se_notch))
  }
  # Compute A and b for S_a = A N^b
  params <- lm(logS ~ poly(logN, 1, raw=TRUE))
  A <- exp(summary(params)$coefficients[1,1])
  b <- summary(params)$coefficients[2,1]
  N_Sar1 <- (Snet/A)^(1/b)
  S_trace <- Snet
  Ncurve_trace <- c(1,N_Sar1,N_Sar1)
  Scurve_trace <- c(Snet,Snet,round_any(Se_notch, 10^floor(log10(Se_notch)), f = floor))

  df <- data.frame(S = c(Scurve_unnotched,Scurve_notched), S1 = c(Scurve_unnotched,rep(NA,4)), S2 = c(rep(NA,4),Scurve_notched), N = c(Ncurve,Ncurve), notchstate = c(rep("unnotched",4),rep("notched",4)))
  df2 <- data.frame(St = Scurve_trace, Nt = Ncurve_trace)

  plotout<-ggplot() +
    geom_path(data=df, aes(N,S, colour = notchstate), size = 0.9) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab(Xlab) +
    ylab(Ylab)
  plotout <- plotout + geom_line(data=df2, aes(Nt,St), colour = "blue", size = 0.9, linetype = "dashed")

  return(list(plot1 = plotout, Kt = Kt, q = q, Kf = Kf, Kfp = Kfp, unnotched_S = list(at1 = S1_notch, at1000 = S1000, at1000000 = Se), notched_S = list(at1 = S1_notch, at1000 = S1000_notch, at1000000 = Se_notch), Snet = Snet, A = A, b = b, Nf = N_Sar1))
}
