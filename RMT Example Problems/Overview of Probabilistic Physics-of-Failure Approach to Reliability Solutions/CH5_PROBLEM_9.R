# CHAPTER 5 PROBLEM 9
# Reuel Smith
# =================================================
# We are interested to develop the probability of detection (POD) and measurement error
# model associated with corrosion pit detection using optical microscopy.

#            No of pits     Poisson(x) ~ Pr(ρ|x)
# Density = ------------ = ----------------------
#            unit area           4 mm x 4mm

# Overall area of 10 square centimeter sample of high carbon steel is immersed
# in a 6% ferric chloride solution with 1% hydrochloric acid (HCl) to stabilize
# the pH at a temperature of 45 °C for 72 hours.

# With a very high precision technique the true pit sizes were recorded on a
# 4mm x 4mm section of the sample after the test. Then inspectors were asked
# to use a 50X optical microscope to detect the same pits on the same 4 mm x 4mm area.
# 34 pits identified

# Objective: Develop POD(y) and measurement error y = m yhat + c models
# Assume prior flaw number is poisson distribution where pit density  ρ (per cm^2) is defined NOR(200,50)
# truncated to  ρ > 0 and where pit size prior is LOGN(log(30), 0.8)

a_lth <- 0.254  # State the threshold

# Enter the full data set including measured, actual size, and censored
# Include non-detected measurements as 'NA' for completeness (pit size in μm)
meas_a <- c(33, 33, 24, 103, 59, NA, 63, 23, 27, 39, NA, NA, 24, NA, 47, 17,
            10, 95, 33, 38, NA, 20, 33, NA, 32, 25, 29, 20, 138, 22, 14, 72)
actual_a <- c(27, 25, 25, 106, 55, 21, 46, 21, 21, 36, 6, 8, 19, 13, 39, 15,
              10, 97, 25, 32, 10, 19, 33, 12, 41, 19, 24, 19, 132, 30, 13, 58)
censored_a <- c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0,
                1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1)

# Bayesian update of density
# lambda = density*A
output.CH5.P9.a <- distribution.BAYESest(200,dist="Poisson",TTF=0.16,Tc=NULL,N=32,Kf=NULL,confid = 0.95,priors = "normal(200,50) T[0 , ]",nsamples=100000,burnin=1000,nchains=4)

output.CH5.P9.b <- distribution.BAYESest(c(3,0.8),dist="Lognormal",TTF=actual_a,confid = 0.95,priors = "lognormal(3.401197,0.8)",nsamples=100000,burnin=1000,nchains=4)
probplot.logn(data = cbind(actual_a,rep(1,32),rep(1,32)), pp = "Blom", xlabel1 = "Size", MLE_i = 1)

log(30)*((0.8^2)/(0.6905921^2 + 0.8^2)) + 3.226388*((0.6905921^2)/(0.6905921^2 + 0.8^2))
1/sqrt((1/0.8^2) + (1/0.6905921^2))

data_POD_P9 <- cbind(meas_a,actual_a,censored_a)
POD.calc(data_POD_P9,dist = "Logistic",alth = 0.01,xlabel1 = "Pit Size (μm)",ylabel1 = "POD",unitlabel = "μm")
# -7.691532
# 19.38306
POD.calc(data_POD_P9,dist = "Loglogistic",alth = 0.01,xlabel1 = "Pit Size (μm)",ylabel1 = "POD",unitlabel = "μm")
# -7.601296
# 19.20259
POD.calc(data_POD_P9,dist = "Lognormal",alth = 0.01,xlabel1 = "Pit Size (μm)",ylabel1 = "POD",unitlabel = "μm")
# -7.550974
# 19.10195

priors.P9 <- c("normal(-16.5,3.75)","normal(4.2,1.5)")
distparam <-"real alpha_1; real alpha_2;"
distpriors<-paste(c("alpha_1 ~ ",priors.P9[1],";","alpha_2 ~ ",priors.P9[2],";"),collapse = "")

loglik <- paste(c("target += bernoulli_lpmf(D | (1./(1 + exp(-(alpha_1 + alpha_2.*log(a))))));"),collapse = "")
params <- paste(c(distparam),collapse = " ")
paramsvec <- c("alpha_1","alpha_2")
priors <- paste(c(distpriors),collapse = " ")

block1 <- "data {int n; vector[n] a; array[n] int D;}"
datablock <- list(n = length(actual_a), a = actual_a, D = censored_a)

block2 <- paste(c("parameters {",params,"}"),collapse = " ")

block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

stanlscode <- paste(c(block1,block2,block3),collapse=" ")

stanlsfile <- write_stan_file(stanlscode)
print(stanlsfile)
# Generate initial list (one list per chain)
pt_est <- c(1,1)
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

# Compute hit-or-miss
Actual_data <- actual_a
True_Flaw_individual <- sort(unique(actual_a)) # get the true size points
Hit_vector <- rep(0,length(True_Flaw_individual))
Miss_vector <- rep(0,length(True_Flaw_individual))
Total_vector <- rep(0,length(True_Flaw_individual))

for(i in 1:length(True_Flaw_individual)){ # Hit and miss counts against real  size
  Hit_vector[i] <- sum(censored_a[which(Actual_data == True_Flaw_individual[i])])
  Miss_vector[i] <- length(which(Actual_data == True_Flaw_individual[i])) - sum(censored_a[which(Actual_data == True_Flaw_individual[i])])
}
Total_vector <- Hit_vector + Miss_vector # Total checks per size
POD_point.est <- Hit_vector/Total_vector # Point estimate for POD
logit <- log(POD_point.est/(1 - POD_point.est)) # Compute log-odds or logit

# Compute Credible Intervals for POD
alpha <- 0.05
percentiles <- c(0, 1) + (c(1, -1) * (alpha / 2))
for(i in 1:length(xplot)){
  POD.CI.0 <- as.vector(quantile((1./(1 + exp(-(extract(fit,c("alpha_1"))$alpha_1 + stats.mean.sd[2,1]*log(xplot[i]))))),
           probs = percentiles))
  if(i == 1){
    POD.CI <- POD.CI.0
  }
  if(i > 1){
    POD.CI <- rbind(POD.CI,POD.CI.0)
  }
}


# Generate plot
xplot <- linspace(0.1,max(Actual_data),1000)
yplot_MLE <- (1./(1 + exp(-(stats.mean.sd[1,1] + stats.mean.sd[2,1]*log(xplot)))))

df_HitMiss <- data.frame(x = c(Actual_data,True_Flaw_individual), y = c(censored_a, POD_point.est), H.or.M = c(rep("Hit or Miss Data",length(actual_a)),rep("POD Point Estimates",length(True_Flaw_individual))))
df_POD_curve <- data.frame(x = c(xplot), y = c(yplot_MLE),POD.curve = c(rep("POD",1000)))
df_POD_CI <- data.frame(x = c(xplot), ymin = as.vector(POD.CI[,1]), ymax = as.vector(POD.CI[,2]),POD.curve = c(rep("POD CI",1000)))

plotout <- ggplot() + geom_point(data = df_HitMiss, aes(x, y, shape = H.or.M), size = 2.5, colour = "black") +
  geom_path(data = df_POD_curve, aes(x, y, linetype = POD.curve), linewidth = 0.4, colour = "black") +
  geom_ribbon(data = df_POD_CI, aes(x=x, ymin = ymin, ymax = ymax, fill = POD.curve ),alpha = 0.25) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab("Pit Size (μm)") +
  ylab("POD")

# Density plot
posterior_alpha_1 <- density(extract(fit,c("alpha_1"))$alpha_1)
posterior_alpha_2 <- density(extract(fit,c("alpha_2"))$alpha_2)
df_posterior <- data.frame(x = c(posterior_alpha_1$x,posterior_alpha_2$x),
                           ymin = rep(0,(length(posterior_alpha_1$x)+length(posterior_alpha_2$x))),
                           ymax = c(posterior_alpha_1$y,posterior_alpha_2$y),
                           distlabel = c(rep("alpha_1",length(posterior_alpha_1$x)),rep("alpha_2",length(posterior_alpha_2$x))))

# Density plot for Exponential rate parameter posterior
plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  facet_wrap(~distlabel, dir="v", scales = "free") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab(" ") +
  ylab("density")

# Measurement Error
MeasError.calc(data_POD_P9,xlabel1 = "True Pit Size (μm)",ylabel1 = "Measured Pit Size (μm)",ME = "Additive")
MeasError.calc(data_POD_P9,ylabel1 = "Model Pit Size (μm)",xlabel1 = "Measured Pit Size (μm)",ME = "Multiplicative")

# priors.P9.ME <- c("lognormal(2.533792,0.07538997)","lognormal(0.01042184,0.03733932)","normal(1.21408,2.032681)")
priors.P9.ME <- c("uniform(0.01,5)","uniform(0.1,2)","uniform(-5,1)")
params <- "real<lower=0> sigma; real m; real c;"
paramsvec <- c("sigma","m","c")
distpriors<-paste(c("sigma ~ ",priors.P9.ME[1],";","m ~ ",priors.P9.ME[2],";","c ~ ",priors.P9.ME[3],";"),collapse = "")

lifeF <- "m*log(De) + c"

loglik <- paste(c("target += lognormal_lpdf(Dtrue |",lifeF,", sigma);"),collapse = "")
priors <- distpriors
outputparamset <- c("\U03C3","m","c")

# BLOCK 1
block1 <- "data {int<lower=0> n; vector[n] Dtrue; vector[n] De;}"
datablock <- list(n = length(actual_a[which(is.na(meas_a)==FALSE)]),
                  Dtrue = actual_a[which(is.na(meas_a)==FALSE)],
                  De = meas_a[which(is.na(meas_a)==FALSE)])

# BLOCK 2
block2 <- paste(c("parameters {",params,"}"),collapse = " ")

# BLOCK 3
block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

stanlscode <- paste(c(block1,block2,block3),collapse=" ")

stanlsfile <- write_stan_file(stanlscode)
print(stanlsfile)

pt_est <- c(0.5,0.98,-0.038)
# Generate initial list (one list per chain)
names(pt_est) <- paramsvec
pt_estlist <- as.list(pt_est)
init_pt_est <- vector("list",4)
for(i in 1:4){
  init_pt_est[[i]] <- pt_estlist
}
# Build or compile Stan code to C++
lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
fit <- sampling(lsmod, data = datablock, iter = 50000, warmup = 1000, init = init_pt_est)
stats.mean.sd <- summary(fit)$summary[,c(1,3)]
stats.Rhat <- rhat(fit)
confidbounds <- mcmc_intervals_data(data.frame(extract(fit, paramsvec)),prob_outer = 0.95)
outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec),1],unname(stats.mean.sd)[1:length(paramsvec),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec)]), nrow = length(paramsvec), ncol = 6, byrow = FALSE,dimnames = list(paramsvec,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))

cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
print(outputtable)
cat(c("\n"),sep = "")

# Construct new ME plot and Bayesian CI
# Lines and Confidence
measured_line <- linspace(0.1,max(actual_a),100)
mean_model_line <- exp(stats.mean.sd[3,1] + stats.mean.sd[2,1]*log(measured_line))
# var_model_line <- (log(measured_line)^2)*varcov1[1,1] + varcov1[2,2]
crit <- qnorm((1 + 0.95)/2)
# Compute Credible Intervals for POD
alpha <- 0.05
percentiles <- c(0, 1) + (c(1, -1) * (alpha / 2))
for(i in 1:length(measured_line)){
  Model.Fit.0 <- as.vector(quantile(exp(extract(fit,c("c"))$c + extract(fit,c("m"))$m*log(measured_line[i])),
                                 probs = percentiles))
  if(i == 1){
    Model.Fit.CI <- Model.Fit.0
  }
  if(i > 1){
    Model.Fit.CI <- rbind(Model.Fit.CI,Model.Fit.0)
  }
}
# lower_model_line <- exp(log(mean_model_line) - crit * sqrt(var_model_line) + qnorm(((1-confid)/2),0,theta.hat1[1]))
# upper_model_line <- exp(log(mean_model_line) + crit * sqrt(var_model_line) + qnorm((1-((1-confid)/2)),0,theta.hat1[1]))

# return(theta.hat1)
df_data <- data.frame(x = meas_a[which(is.na(meas_a)==FALSE)], y = actual_a[which(is.na(meas_a)==FALSE)], Data.points = rep("Data",length(which(is.na(meas_a)==FALSE))))
df_line <- data.frame(x = c(measured_line,measured_line,measured_line), y = c(mean_model_line,Model.Fit.CI[,1],Model.Fit.CI[,2]), best.fit = c(rep("Multiplicative ME Model",100),rep("Lower Bayesian CI",100),rep("Upper Bayesian CI",100)))
# df_line <- data.frame(x = c(measured_line), y = c(mean_model_line), best.fit = c(rep("Multiplicative ME Model",100)))

plotout <- ggplot() + geom_path(data = df_line, aes(x, y, linetype = best.fit), linewidth = 0.4, colour = "black") +
  geom_point(data = df_data, aes(x, y, shape = Data.points), colour = "black") +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),size = .4)) +
  scale_x_continuous(expand=c(0, 0), limits = c(0, max(actual_a))) +
  scale_y_continuous(expand=c(0, 0), limits = c(0, max(actual_a))) +
  xlab("Measured Pit Size (μm)") +
  ylab("Model Pit Size (μm)") +
  scale_linetype_discrete("Best Fit") +
  scale_shape_discrete("Data Points")
plotout

# Density plot
posterior_sigma <- density(extract(fit,c("sigma"))$sigma)
posterior_m <- density(extract(fit,c("m"))$m)
posterior_c <- density(extract(fit,c("c"))$c)
df_posterior <- data.frame(x = c(posterior_sigma$x,posterior_m$x,posterior_c$x),
                           ymin = rep(0,(length(posterior_sigma$x)+length(posterior_m$x)+length(posterior_c$x))),
                           ymax = c(posterior_sigma$y,posterior_m$y,posterior_c$y),
                           distlabel = c(rep("sigma",length(posterior_sigma$x)),rep("m",length(posterior_m$x)),rep("c",length(posterior_c$x))))

# Density plot for Exponential rate parameter posterior
plot3_density.ME <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  facet_wrap(~distlabel, dir="v", scales = "free") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab(" ") +
  ylab("density")
plot3_density.ME

