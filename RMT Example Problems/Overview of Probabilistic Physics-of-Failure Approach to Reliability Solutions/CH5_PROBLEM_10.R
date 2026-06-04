# CHAPTER 5 PROBLEM 10
# Reuel Smith
# =================================================
# STEP 1: Set up data and posterior samples
# Posterior samples for POD model
POST_SAMPLE_POD <- data.frame(extract(fit.1, c("alpha_1","alpha_2")))

# Posterior samples for measurement error model
POST_SAMPLE_ME <- data.frame(extract(fit, c("m","c","sigma")))

# Measured Time and Measured pit radius data
Time.P10 <- c(rep(c(48,72,96),25),
              rep(c(72,96),21))
Meas_R.P10 <- c(4, 16, 31,
                8, 18, 38,
                10, 20, 42,
                11, 21, 47,
                11, 22, 49,
                12, 23, 52,
                13, 25, 54,
                13, 25, 55,
                13, 25, 55,
                14, 26, 56,
                15, 26, 60,
                16, 26, 61,
                16, 27, 61,
                18, 27, 62,
                18, 28, 62,
                19, 29, 64,
                19, 29, 66,
                19, 29, 66,
                20, 30, 68,
                20, 31, 69,
                20, 33, 69,
                21, 33, 70,
                22, 34, 72,
                23, 35, 73,
                25, 35, 81,
                36, 81,
                37, 83,
                37, 85,
                38, 87,
                38, 89,
                39, 89,
                39, 104,
                40, 105,
                41, 112,
                42, 114,
                45, 115,
                46, 120,
                46, 134,
                48, 136,
                48, 137,
                48, 137,
                49, 181,
                49, 187,
                50, 197,
                54, 204,
                55, 280)
Unit.P10 <- c(rep("Unit 1",3),rep("Unit 2",3),rep("Unit 3",3),rep("Unit 4",3),rep("Unit 5",3),
              rep("Unit 6",3),rep("Unit 7",3),rep("Unit 8",3),rep("Unit 9",3),rep("Unit 10",3),
              rep("Unit 11",3),rep("Unit 12",3),rep("Unit 13",3),rep("Unit 14",3),rep("Unit 15",3),
              rep("Unit 16",3),rep("Unit 17",3),rep("Unit 18",3),rep("Unit 19",3),rep("Unit 20",3),
              rep("Unit 21",3),rep("Unit 22",3),rep("Unit 23",3),rep("Unit 24",3),rep("Unit 25",3),
              rep("Unit 26",2),rep("Unit 27",2),rep("Unit 28",2),rep("Unit 29",2),rep("Unit 30",2),
              rep("Unit 31",2),rep("Unit 32",2),rep("Unit 33",2),rep("Unit 34",2),rep("Unit 35",2),
              rep("Unit 36",2),rep("Unit 37",2),rep("Unit 38",2),rep("Unit 39",2),rep("Unit 40",2),
              rep("Unit 41",2),rep("Unit 42",2),rep("Unit 43",2),rep("Unit 44",2),rep("Unit 45",2),
              rep("Unit 46",2))
# STEP 2: Randomly sample the true pit size based on the posterior measurement error model and given data.
True_R.P10 <- rep(0,length(Meas_R.P10))   # Sampled from posterior and data

N <- length(POST_SAMPLE_ME$m)

for(i in 1:length(Meas_R.P10)){
  # Pull measurement error parameters randomly but keep them within the same unit
  if(i == 1){ # For first unit
    i2 <- round(runif(1,1,N)) # Sample from ME fit
    True_R.P10[i] <- rlnorm(1,POST_SAMPLE_ME$m[i2]*log(Meas_R.P10[i]) + POST_SAMPLE_ME$c[i2],POST_SAMPLE_ME$sigma[i2])
  }
  if(i > 1 && isTRUE(Unit.P10[i-1] != Unit.P10[i])==TRUE){ # For changed unit number
    i2 <- round(runif(1,1,N)) # Sample from ME fit
    True_R.P10[i] <- rlnorm(1,POST_SAMPLE_ME$m[i2]*log(Meas_R.P10[i]) + POST_SAMPLE_ME$c[i2],POST_SAMPLE_ME$sigma[i2])
  }
  if(i > 1 && isTRUE(Unit.P10[i-1] != Unit.P10[i])==FALSE){ # For same unit number
    plow <- plnorm(True_R.P10[i-1],POST_SAMPLE_ME$m[i2]*log(Meas_R.P10[i]) + POST_SAMPLE_ME$c[i2],POST_SAMPLE_ME$sigma[i2])
    True_R.P10[i] <- qlnorm(runif(1,plow,1),POST_SAMPLE_ME$m[i2]*log(Meas_R.P10[i]) + POST_SAMPLE_ME$c[i2],POST_SAMPLE_ME$sigma[i2])
  }
}

# STEP 3: Determine the detectability of these samples based on the posterior POD and a random detection probability draw per sample.
Detected.P10 <- rep(1,length(True_R.P10))   # Default detected
N2 <- length(POST_SAMPLE_POD$alpha_1)

for(i in 1:length(Meas_R.P10)){
  # Pull POD parameters randomly but again keep them within the same unit
  if(i == 1){ # For first unit
    i2 <- round(runif(1,1,N2)) # Sample from POD fit
    alp_1.0 <- POST_SAMPLE_POD$alpha_1[i2]
    alp_2.0 <- POST_SAMPLE_POD$alpha_2[i2]
  }
  if(i > 1 && isTRUE(Unit.P10[i-1] != Unit.P10[i])==TRUE){ # For changed unit number
    i2 <- round(runif(1,1,N2)) # Sample from POD fit
    alp_1.0 <- POST_SAMPLE_POD$alpha_1[i2]
    alp_2.0 <- POST_SAMPLE_POD$alpha_2[i2]
  }

  POD.val <- 1/(1 + exp(-(alp_1.0 + alp_2.0*True_R.P10[i])))
  if(runif(1) > POD.val){
    Detected.P10[i] <- NA
  }
}

# STEP 4: Use this simulated true data and given time to run Bayesian analysis for   assuming a non-informative prior.
priors.P10 <- c("normal(60000000,20000000)","uniform(0.1,100)")
distparam <-"real DH; real<lower=0> sigma_t; "
distpriors<-paste(c("DH ~ ",priors.P10[1],";","sigma_t ~ ",priors.P10[2],";"),collapse = "")

loglik <- paste(c("target += lognormal_lpdf(Size |0.5.*log(4095.871.*Time*exp(-DH/2645099)), sigma_t);"),collapse = "")
params <- paste(c(distparam),collapse = " ")
paramsvec <- c("DH","sigma_t")
priors <- paste(c(distpriors),collapse = " ")

block1 <- "data {int<lower=0> n; vector[n] Time; vector[n] Size;}"

datablock <- list(n = length(Detected.P10), Time = Time.P10, Size = True_R.P10)

block2 <- paste(c("parameters {",params,"}"),collapse = " ")

block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

stanlscode <- paste(c(block1,block2,block3),collapse=" ")

stanlsfile <- write_stan_file(stanlscode)
print(stanlsfile)
# Generate initial list (one list per chain)
pt_est <- c(60000000,0.5)
conf.level <- 0.95
names(pt_est) <- paramsvec
pt_estlist <- as.list(pt_est)
init_pt_est <- vector("list",4)
for(i in 1:4){
  init_pt_est[[i]] <- pt_estlist
}

conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))

lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
fit.1 <- sampling(lsmod, data = datablock, iter = 50000, warmup = 1000, init = pt_est)
stats.mean.sd <- summary(fit.1)$summary[,c(1,3)]
stats.Rhat <- rhat(fit.1)
confidbounds <- mcmc_intervals_data(data.frame(extract(fit.1, paramsvec)),prob_outer = 0.95)
outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec),1],unname(stats.mean.sd)[1:length(paramsvec),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec)]), nrow = length(paramsvec), ncol = 6, byrow = FALSE,dimnames = list(paramsvec,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))

cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
print(outputtable)
cat(c("\n"),sep = "")



# Compute Credible Intervals for Degradation Model r(t)
alpha <- 0.05
percentiles <- c(0, 1) + (c(1, -1) * (alpha / 2))
for(i in 1:length(xplot)){
  R.CI.0 <- as.vector(quantile(sqrt(4095.871*xplot[i]*exp(-extract(fit.1,c("DH"))$DH/2645099)),
                               probs = percentiles))
  if(i == 1){
    R.CI <- R.CI.0
  }
  if(i > 1){
    R.CI <- rbind(R.CI,R.CI.0)
  }
}


# Generate plot
xplot <- linspace(0.1,max(True_R.P10),1000)
yplot_MLE <- sqrt(4095.871*xplot*exp(-stats.mean.sd[1,1]/2645099))

df_data <- data.frame(x = c(Time.P10), y = c(True_R.P10), Data.set = rep("Simulated True Data",length(Time.P10)))
df_R_curve <- data.frame(x = c(xplot), y = c(yplot_MLE),Pit.curve = c(rep("Pit Size Model Mean",1000)))
df_R_CI <- data.frame(x = c(xplot), ymin = as.vector(R.CI[,1]), ymax = as.vector(R.CI[,2]),Pit.curve = c(rep("Pit Size CI",1000)))

plotout.3 <- ggplot() + geom_point(data = df_data, aes(x, y, shape = Data.set), size = 2.5, colour = "black") +
  geom_path(data = df_R_curve, aes(x, y, linetype = Pit.curve), linewidth = 0.4, colour = "black") +
  geom_ribbon(data = df_R_CI, aes(x=x, ymin = ymin, ymax = ymax, fill = Pit.curve ),alpha = 0.25) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab("Time (Hours)") +
  ylab("Pit Size (μm)")
plotout.3

# Density plot
posterior_DH <- density(extract(fit.1,c("DH"))$DH)
posterior_sigma_t <- density(extract(fit.1,c("sigma_t"))$sigma_t)
df_posterior <- data.frame(x = c(posterior_DH$x,posterior_sigma_t$x),
                           ymin = rep(0,(length(posterior_DH$x)+length(posterior_sigma_t$x))),
                           ymax = c(posterior_DH$y,posterior_sigma_t$y),
                           distlabel = c(rep("Activation Energy, ΔH (J/mol)",length(posterior_DH$x)),rep("σ_t",length(posterior_sigma_t$x))))

# Density plot for Exponential rate parameter posterior
plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  facet_wrap(~distlabel, dir="v", scales = "free") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab(" ") +
  ylab("density")
plot3_density
