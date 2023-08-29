# Stress-Strain Parameters
# Developed by Reuel Smith, 2022

stress_strain.params <- function(dat,E,stressunits,options){
  # dat is entered as a list made up of stress, strain, and cycles in that order
  library(pracma)
  library(nls.multstart)
  library(ggplot2)

  if(missing(dat)){
    stop('Enter data for stress, strain, and cycles (in list form).')
  }
  if(is.list(dat)==FALSE){
    stop('Enter data as list(stress,strain,cycles).')
  }
  if(missing(E)){
    stop('Enter modulus of elasticity E (in MPa or ksi).')
  }

  # Check units and set up axis labels
  if(missing(stressunits) || stressunits == 1){
    stressunits <- c("MPa")
  }
  if(stressunits == 2){
    stressunits <- c("ksi")
  }

  # Separate the data if in list form
  stress <- dat[[1]]
  strain <- dat[[2]]
  cycles <- dat[[3]]

  # Check units and set up axis labels
  Xlab1 <- paste(c("Total Strain, (eps_tot)"),collapse="", sep = "_")
  Ylab1 <- paste(c("Alternating Stress, S_a (",stressunits,")"),collapse="", sep = "_")
  Xlab2 <- paste(c("Reversals to Failure, 2N_f"),collapse="", sep = "_")
  Ylab2 <- paste(c("Strain amplitude, eps_a"),collapse="", sep = "_")

  # Refer to options.  If no other relationship is stated,
  # default to Coffin-Mason Relationship for evaluation
  # ========================================
  # Elastic Strain and sig_f' and b
  # ========================================
  paramset1 <- lm(log(stress) ~ poly(log(2*cycles), 1, raw=TRUE))
  sig_f <- exp(summary(paramset1)$coefficients[1,1])
  b <- summary(paramset1)$coefficients[2,1]
  # ========================================
  # Plastic Strain and eps_f', c, K, and n
  # =======================================
  strain_p <- strain - (stress/E)
  paramset2 <- lm(log(strain_p) ~ poly(log(2*cycles), 1, raw=TRUE))
  eps_f <- exp(summary(paramset2)$coefficients[1,1])
  c <- summary(paramset2)$coefficients[2,1]

  paramset3 <- lm(log(strain_p) ~ poly(log(stress), 1, raw=TRUE))
  n <- 1/summary(paramset3)$coefficients[2,1]
  K <- exp(-n*summary(paramset3)$coefficients[1,1])
  # =======================================
  # Modifier Options
  # =======================================
  if(isFALSE(missing(options))){
    if(length(options$mean_stress_corr) == 1 && length(options$Sm) == 0){
      stop('Please enter a mean stress Sm.')
    }
    if(length(options$mean_stress_corr) == 1 && length(options$Sm) == 1){
      Sm <- options$Sm
    }
  }
  # =======================================
  # Plotting Output
  # =======================================
  # For stress-strain curve
  stressline <- linspace(0,max(stress),100)
  strainline <- (stressline/E) + (stressline/K)^(1/n)
  df1 <- data.frame(strain_ = c(strainline,rep(NA,length(strain))), stress_ = c(stressline,rep(NA,length(stress))), strain_data = c(rep(NA,length(strainline)),strain), stress_data = c(rep(NA,length(stressline)),stress), data_points = c(rep("curve",length(stressline)),rep("data",length(stress))))

  plotout1<-ggplot() +
    geom_point(data=df1, aes(strain_data,stress_data), colour = 'red', size = 1.9) +
    geom_line(data=df1, aes(strain_,stress_), colour = "black", size = 0.9) +
    xlab(Xlab1) +
    ylab(Ylab1)

  # For Strain-Life curve
  revline <- exp(linspace(log(1),log(10^(floor(log10(2*max(cycles)))+2)),100))
  elasticstrain <- (sig_f/E)*(revline^b)
  plasticstrain <- eps_f*(revline^c)
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && (options$mean_stress_corr == "Morrow" || options$mean_stress_corr == "ModifiedMorrow") ){
    elasticstrain <- ((sig_f - Sm)/E)*(revline^b)
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "ModifiedMorrow" ){
    plasticstrain <- eps_f*(((sig_f - Sm)/sig_f)^(c/b))*(revline^c)
  }
  # =============================================================
  # EVALUATION OF SWT COMES AFTER HYSTERESIS LOOP IF IT EXISTS
  # =============================================================
  # Provide Hysteresis Plot if loadconditions optional input exists
  strain_amp <- logical(0)
  stress_amp <- logical(0)
  if(isFALSE(missing(options))){
    if(length(options$stressrange)>0){
      hysteresisoutput <- hysteresisloop.plot(E,K,n,stressunits,list(stressrange = options$stressrange))
    }
    if(length(options$maxstress)>0){
      hysteresisoutput <- hysteresisloop.plot(E,K,n,stressunits,list(maxstress = options$maxstress))
    }
    if(length(options$minstress)>0){
      hysteresisoutput <- hysteresisloop.plot(E,K,n,stressunits,list(minstress = options$minstress))
    }
    if(length(options$strainrange)>0){
      hysteresisoutput <- hysteresisloop.plot(E,K,n,stressunits,list(strainrange = options$strainrange))
    }
    if(length(options$maxstrain)>0){
      hysteresisoutput <- hysteresisloop.plot(E,K,n,stressunits,list(maxstrain = options$maxstrain))
    }
    if(length(options$minstrain)>0){
      hysteresisoutput <- hysteresisloop.plot(E,K,n,stressunits,list(minstrain = options$minstrain))
    }
    strain_amp <- (hysteresisoutput$maxstrain - hysteresisoutput$minstrain)*0.5
    stress_amp <- (hysteresisoutput$maxstress - hysteresisoutput$minstress)*0.5
    plotout3 <- hysteresisoutput[[1]]
  }
  # =============================================================
  # EVALUATE SWT NOW
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT" ){
    # SWT Mean Stress Correction
    # First check and see if there is a stress amplitude
    if(length(stress_amp)==0){
      stop('SWT can only be computed with current load conditions entered in options.')
    }
    stressmax_corr <- stress_amp + Sm
    elasticstrain <- ((sig_f^2)/(E*stressmax_corr))*(revline^(2*b))
    plasticstrain <- ((eps_f*sig_f)/stressmax_corr)*(revline^(b+c))
  }
  # EVALUATE Walker NOW
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Walker" && length(options$gam) == 1){
    # Walker Mean Stress Correction
    # First check and see if there is a stress amplitude
    if(length(stress_amp)==0){
      stop('Walker can only be computed with current load conditions entered in options.')
    }
    # Calculate load ratio
    R <- hysteresisoutput$minstress/hysteresisoutput$maxstress
    gam <- options$gam
    elasticstrain <- (sig_f/E)*((0.5*(1 - R))^(1 - gam))*(revline^b)
    plasticstrain <- eps_f*((0.5*(1 - R))^((c/b)*(1 - gam)))*(revline^c)
  }
  totalstrain <- elasticstrain + plasticstrain

  df2 <- data.frame(reversals_dat = c(2*cycles,rep(NA,300)), totstrain_dat = c(strain,rep(NA,300)), reversals = c(rep(NA,length(cycles)),revline,revline,revline), totstrain = c(rep(NA,length(strain)),totalstrain,elasticstrain,plasticstrain), strain_lines = c(rep(NA,length(strain)),rep("Total Strain",length(totalstrain)),rep("Elastic Strain",length(elasticstrain)),rep("Plastic Strain",length(plasticstrain))))
  plotout2<-ggplot() +
    geom_point(data=df2, aes(reversals_dat,totstrain_dat), colour = 'red', size = 1.9) +
    geom_line(data=df2, aes(reversals, totstrain, colour = strain_lines), size = 0.9) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    annotation_logticks() +
    xlab(Xlab2) +
    ylab(Ylab2)

  # =======================================
  # Output Field
  # =======================================
  if(missing(options)==TRUE || (isFALSE(missing(options)) && length(options$mean_stress_corr)==0)){
    # Coffin-Mason
    if(length(strain_amp)==1){
      straintrace <- function(reversal) strain_amp - (sig_f/E)*((reversal)^b) - eps_f*((reversal)^c)
    }
    cat(c("Coffin-Mason Model: \U03B5_tot(N_f) = (\U03C3_f'/E) (2 N_f)^b + \U03B5_f' (2 N_f)^c \n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Morrow" ){
    # Morrow Mean Stress Correction
    if(length(strain_amp)==1){
      straintrace <- function(reversal) strain_amp - ((sig_f - Sm)/E)*((reversal)^b) - eps_f*((reversal)^c)
    }
    cat(c("Morrow Mean Stress Correction Model: \U03B5_tot(N_f) = ((\U03C3_f' - \U03C3_m)/E) (2 N_f)^b + \U03B5_f' (2 N_f)^c \n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "ModifiedMorrow" ){
    # Modified Morrow Mean Stress Correction
    if(length(strain_amp)==1){
      straintrace <- function(reversal) strain_amp - ((sig_f - Sm)/E)*((reversal)^b) - eps_f*(((sig_f - Sm)/sig_f)^(c/b))*((reversal)^c)
    }
    cat(c("Modified Morrow Mean Stress Correction Model: \U03B5_tot(N_f) = ((\U03C3_f' - \U03C3_a)/E) (2 N_f)^b + \U03B5_f' ((\U03C3_f' - \U03C3_a)/\U03C3_f')^(c/b) (2 N_f)^c \n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT" ){
    # SWT Mean Stress Correction
    if(length(strain_amp)==1){
      straintrace <- function(reversal) strain_amp - (1/stressmax_corr)*((sig_f^2)/E)*((reversal)^(2*b)) - (1/stressmax_corr)*eps_f*sig_f*((reversal)^(b+c))
    }
    cat(c("SWT Mean Stress Correction Model: \U03C3_max \U03B5_tot(N_f) = ((\U03C3_f'^2)/E) (2 N_f)^2b + \U03C3_f' \U03B5_f' (2 N_f)^(b+c) \n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Walker" ){
    # Walker Mean Stress Correction
    if(length(strain_amp)==1){
      straintrace <- function(reversal) strain_amp - (sig_f/E)*((0.5*(1 - R))^(1 - gam))*((reversal)^b) - eps_f*((0.5*(1 - R))^((c/b)*(1 - gam)))*((reversal)^c)
    }
    cat(c("Walker Mean Stress Correction Model: \U03B5_tot(N_f) = (\U03C3_f'/E) ((1 - R)/2)^(1 - \U03B3) (2 N_f)^b + \U03B5_f' ((1 - R)/2)^((c(1 - \U03B3))/b) (2 N_f)^c \n\n"),sep = "")
  }
  if(length(strain_amp)==1){
    rev_2N_trace <- findzeros(straintrace, 10, max(revline))
    Life_trace <- rev_2N_trace/2
    rev_2Ncurve_trace <- c(min(revline),rev_2N_trace,rev_2N_trace)
    straincurve_trace <- c(strain_amp,strain_amp,min(plasticstrain))
    df3 <- data.frame(reversals_trace = rev_2Ncurve_trace, strain_trace = straincurve_trace)
    plotout2 <- plotout2 + geom_line(data=df3, aes(reversals_trace,strain_trace), colour = "blue", size = 0.9, linetype = "dashed")
  }

  cat(c("Total strain model \U03B5_tot(\U03C3) = \U03C3_a/",E," + (\U03C3_a/",signif(K, digits=6),")^1/",signif(n, digits=4)," or\n"),sep = "")
  if(missing(options)==TRUE || (isFALSE(missing(options)) && length(options$mean_stress_corr)==0)){
    # Coffin-Mason
    cat(c("                   \U03B5_tot(N_f) = ",signif(sig_f/E, digits=4),"(2 N_f)^",signif(b, digits=4)," + ",signif(eps_f, digits=4),"(2 N_f)^",signif(c, digits=4),"\n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Morrow" ){
    # Morrow Mean Stress Correction
    cat(c("                   \U03B5_tot(N_f) = ",signif((sig_f-Sm)/E, digits=4),"(2 N_f)^",signif(b, digits=4)," + ",signif(eps_f, digits=4),"(2 N_f)^",signif(c, digits=4),"\n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "ModifiedMorrow" ){
    # Modified Morrow Mean Stress Correction
    cat(c("                   \U03B5_tot(N_f) = ",signif((sig_f-Sm)/E, digits=4),"(2 N_f)^",signif(b, digits=4)," + ",signif(eps_f*(((sig_f-Sm)/sig_f)^(c/b)), digits=4),"(2 N_f)^",signif(c, digits=4),"\n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "SWT" ){
    # SWT Mean Stress Correction
    cat(c("                   \U03B5_tot(N_f) = ",signif((sig_f^2)/(E*stressmax_corr), digits=4),"(2 N_f)^",signif(2*b, digits=4)," + ",signif((eps_f*sig_f)/stressmax_corr, digits=4),"(2 N_f)^",signif(b+c, digits=4),"\n\n"),sep = "")
  }
  if(isFALSE(missing(options)) && length(options$mean_stress_corr) == 1 && options$mean_stress_corr == "Walker"){
    # Walker Mean Stress Correction
    if(R == -1){
      cat(c("                   \U03B5_tot(N_f) = ",signif(sig_f/E, digits=4)," "," (2 N_f)^",signif(b, digits=4)," + ",signif(eps_f, digits=4)," (2 N_f)^",signif(c, digits=4),"\n\n"),sep = "")
    } else{
      cat(c("                   \U03B5_tot(N_f) = ",signif(sig_f/E, digits=4)," ",signif(0.5*(1 - R), digits=4),"^",signif(1 - gam, digits=4)," (2 N_f)^",signif(b, digits=4)," + ",signif(eps_f, digits=4)," ",signif(0.5*(1 - R), digits=4),"^",signif((c/b)*(1 - gam), digits=4)," (2 N_f)^",signif(c, digits=4),"\n\n"),sep = "")
    }
  }

  cat(c("The estimate for cyclic strain coefficient K' = ",K," ",stressunits,".\n"),sep = "")
  cat(c("The estimate for strain hardening exponent n' = ",n,".\n"),sep = "")
  cat(c("The estimate for fatigue strength coefficient \U03C3_f' = ",sig_f," ",stressunits,".\n"),sep = "")
  cat(c("The estimate for strain hardening exponent \U03B5_f' = ",eps_f,".\n"),sep = "")
  cat(c("The estimate for elastic strain exponent b' = ",b,".\n"),sep = "")
  cat(c("The estimate for plastic strain exponent c' = ",c,".\n\n"),sep = "")

  if(missing(options) || (isFALSE(missing(options)) && sum(length(options$stressrange),length(options$maxstress),length(options$minstress),length(options$strainrange),length(options$maxstrain),length(options$minstrain)) == 0)){
    return(list(plot1 = plotout1, plot2 = plotout2, K = K,n = n,sigma_f = sig_f,epsilon_f = eps_f,b = b,c = c))
  }
  if(sum(length(options$stressrange),length(options$maxstress),length(options$minstress),length(options$strainrange),length(options$maxstrain),length(options$minstrain)) > 0){
    cat(c("The estimate for material failure at strain_a = ",strain_amp," is N = ",Life_trace," cycles or 2N_f = ",rev_2N_trace," reversals.\n\n"),sep = "")
    return(list(plot1 = plotout1, plot2 = plotout2, plot3 = plotout3, K = K,n = n,sigma_f = sig_f,epsilon_f = eps_f,b = b,c = c,maxstrain = hysteresisoutput$maxstrain,minstrain = hysteresisoutput$minstrain,maxstress = hysteresisoutput$maxstress,minstress = hysteresisoutput$minstress, life = Life_trace))
  }
}
