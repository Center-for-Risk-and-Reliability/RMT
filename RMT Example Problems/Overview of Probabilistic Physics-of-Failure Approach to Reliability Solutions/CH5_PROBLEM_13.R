# CHAPTER 5 PROBLEM 13
# Reuel Smith
# =================================================
library(reliabilityRMT)
library(nls.multstart)

Cycles.P13 <- rep(c(30000,60000,90000,120000,150000,180000,210000,240000,270000,300000),6)
Crack.Length.P13 <- c(1.6,	2.1,	2.7,	3.8,	5.5,	7.1,	8.4,	11.7,	15.0,	16.5,
                      1.5,	1.9,	2.8	,3.1,	4.6,	6.2,	8.7,	10.7,	14.1,	20.7,
                      1.1,	1.5,	2.1,	3.3,	4.9,	6.6,	8.5,	11.1,	16.1,	21.6,
                      1.2,	1.7,	2.3,	3.4,	5,	6.7,	10.9,	14.4,	23.3,	26.4,
                      1.8,	3.3,	4.5,	5.8,	6.4,	7.15,	8.4,	8.5,	9.8,	14.2,
                      1.2,	2.5,	5.4,	6.1,	7.59,	8.22,	12.1,	15.3,	15.7,	21.4)/1000

# Bayesian Update and Credible Intervals setup
# Use for Jeffrey's Prior
priors.P13 <- c("target += -2 * log(sigma);")
params <- "real<lower=0> sigma; real<lower=0> C;"
paramsvec <- c("sigma","C")
# distpriors<-paste(c("sigma ~ ",priors.P13[1],";","C ~ ",priors.P13[2],";"),collapse = "")
distpriors<-priors.P13

lifeF <- "log(a_0) + N*3.141593*C*(DS^2)"

loglik <- paste(c("target += lognormal_lpdf(a |",lifeF,", sigma);"),collapse = "")
priors <- distpriors
outputparamset <- c("\U03C3","C")

# BLOCK 1
block1 <- "data {int<lower=0> n; vector[n] N; vector[n] a; real DS; real a_0;}"
datablock <- list(n = length(Cycles.P13),
                  N = Cycles.P13,
                  a = Crack.Length.P13,
                  DS = 200,
                  a_0 = 0.001)

# BLOCK 2
block2 <- paste(c("parameters {",params,"}"),collapse = " ")

# BLOCK 3
block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

stanlscode <- paste(c(block1,block2,block3),collapse=" ")

stanlsfile <- write_stan_file(stanlscode)
print(stanlsfile)

pt_est <- c(0.000000001,0.000000001)
# Generate initial list (one list per chain)
names(pt_est) <- paramsvec
pt_estlist <- as.list(pt_est)
init_pt_est <- vector("list",4)
for(i in 1:4){
  init_pt_est[[i]] <- pt_estlist
}
conf.level <- 0.95
conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))

# Build or compile Stan code to C++
lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
fit <- sampling(lsmod, data = datablock, iter = 100000, warmup = 1000, init = init_pt_est)
stats.mean.sd <- summary(fit)$summary[,c(1,2)]
stats.Rhat <- rhat(fit)
confidbounds <- mcmc_intervals_data(data.frame(extract(fit, paramsvec)),prob_outer = 0.95)
outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec),1],unname(stats.mean.sd)[1:length(paramsvec),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec)]), nrow = length(paramsvec), ncol = 6, byrow = FALSE,dimnames = list(paramsvec,c("Mean","Standard Deviation",conflim_txt[1],"Median",conflim_txt[2],"R\U005E")))

cat(c("Posterior estimates for Bayesian Analysis.\n\n"),sep = "")
print(outputtable)
cat(c("\n"),sep = "")

# Construct new ME plot and Bayesian CI
# Lines and Confidence
Cycles_Plot <- linspace(1,max(Cycles.P13),10000)
mean_model_line <- 0.001*exp(Cycles_Plot*pi*200^2*stats.mean.sd[2,1])
# var_model_line <- (log(measured_line)^2)*varcov1[1,1] + varcov1[2,2]
crit <- qnorm((1 + 0.95)/2)
# Compute Credible Intervals for POD
alpha <- 0.05
percentiles <- c(0, 1) + (c(1, -1) * (alpha / 2))
for(i in 1:length(Cycles_Plot)){
  Model.Fit.0 <- as.vector(quantile(0.001*exp(Cycles_Plot[i]*pi*(200^2)*extract(fit,c("C"))$C),
                                    probs = percentiles))
  if(i == 1){
    Model.Fit.CI <- Model.Fit.0
  }
  if(i > 1){
    Model.Fit.CI <- rbind(Model.Fit.CI,Model.Fit.0)
  }
}

# return(theta.hat1)
df_data <- data.frame(x = Cycles.P13, y = Crack.Length.P13, Data.points = rep("Data",length(Crack.Length.P13)))
df_line <- data.frame(x = c(Cycles_Plot,Cycles_Plot,Cycles_Plot), y = c(mean_model_line,Model.Fit.CI[,1],Model.Fit.CI[,2]), best.fit = c(rep("Simplified Paris Model",10000),rep("Lower Bayesian CI",10000),rep("Upper Bayesian CI",10000)))
# df_line <- data.frame(x = c(measured_line), y = c(mean_model_line), best.fit = c(rep("Multiplicative ME Model",100)))

plotout <- ggplot() + geom_path(data = df_line, aes(x, y, linetype = best.fit), linewidth = 0.4, colour = "black") +
  geom_point(data = df_data, aes(x, y, shape = Data.points), colour = "black") +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),size = .4)) +
  scale_x_continuous(expand=c(0, 0), limits = c(0, max(Cycles.P13))) +
  xlab("Fatigue Cycles") +
  ylab("Crack Length (m)") +
  scale_linetype_discrete("Best Fit") +
  scale_shape_discrete("Data Points")
plotout

# Density plot
posterior_sigma <- density(extract(fit,c("sigma"))$sigma)
posterior_c <- density(extract(fit,c("C"))$C)
df_posterior <- data.frame(x = c(posterior_sigma$x,posterior_c$x),
                           ymin = rep(0,(length(posterior_sigma$x)+length(posterior_c$x))),
                           ymax = c(posterior_sigma$y,posterior_c$y),
                           distlabel = c(rep("sigma",length(posterior_sigma$x)),rep("C",length(posterior_c$x))))

# Density plot for Exponential rate parameter posterior
plot3_density <- ggplot() + geom_ribbon(data = df_posterior, aes(x=x, ymin = ymin, ymax = ymax), fill = "red" ,alpha = 0.5) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  facet_wrap(~distlabel, dir="v", scales = "free") +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab(" ") +
  ylab("density")
plot3_density

# Part b
Df <- 0.0346
# Mean value MCTF
mean_CTF <- (log(Df) - log(0.001))/(pi*200^2*stats.mean.sd[2,1])
# var_model_line <- (log(measured_line)^2)*varcov1[1,1] + varcov1[2,2]
crit <- qnorm((1 + 0.95)/2)
# Compute Credible Intervals for POD
alpha <- 0.05
percentiles <- c(0, 1) + (c(1, -1) * (alpha / 2))
CI_CTF <-  as.vector(quantile((log(Df) - log(0.001))/(pi*200^2*extract(fit,c("C"))$C),
                              probs = percentiles))

for(i in 1:length(Cycles_Plot)){
  Model.Fit.0 <- as.vector(quantile((log(Df) - log(0.001))/(pi*200^2*extract(fit,c("C"))$C),
                                    probs = percentiles))
  Model.Fit.0 <- as.vector(quantile(0.001*exp(Cycles_Plot[i]*pi*(200^2)*extract(fit,c("C"))$C),
                                    probs = percentiles))
  if(i == 1){
    Model.Fit.CI <- Model.Fit.0
  }
  if(i > 1){
    Model.Fit.CI <- rbind(Model.Fit.CI,Model.Fit.0)
  }
}
LSQ.CI_CTF <- qlnorm(percentiles,12.813,0.174)
MLE.CI_CTF <- qlnorm(percentiles,12.981,0.438)

# Box and Whisker Plot
MCTF.box.plot <- data.frame(METHOD = c(rep("Bayesian",length(extract(fit,c("C"))$C)),rep("MLE",1000),rep("LSQ",1000)),
                            value = c((log(Df) - log(0.001))/(pi*200^2*extract(fit,c("C"))$C),rlnorm(1000,12.981,0.438),rlnorm(1000,12.813,0.174)))

ggplot(data = MCTF.box.plot, aes(x=factor(METHOD), y=value, fill=METHOD)) +
  geom_boxplot(outliers = FALSE)+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(axis.title = element_text(family = "serif", size = 10),axis.text = element_text(family = "serif", size = 8), legend.position = c(.06, .74), legend.title = element_text(family = "serif"), legend.text = element_text(family = "serif", size = 8), panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey50")) +
  scale_y_continuous(trans = 'log10') +
  labs(x="Methodololgy",y="Mean Cycle to Failure")
