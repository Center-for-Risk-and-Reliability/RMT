# CHAPTER 5 PROBLEM 12
# Reuel Smith
# =================================================
a_lth <- 0.254  # State the threshold

# Enter the full data set including measured, actual size, and censored
# Include non-detected measurements as 'NA' for completeness
meas_a <- c(2.54, 81.28, 5.08, 99.06, 22.86, NA, 12.7, NA, NA, 12.7, NA, NA, NA, NA, 12.7, 6.35, 48.26, 29.21, 73.66, 7.62,
            109.22, 5.08, NA, NA, 5.08, 25.4, NA, 82.55, 12.7, 101.6, 38.1, 6.35, 76.2, 88.9, 10.16, NA, NA, NA, 27.94, 76.2,
            38.1, 25.4, 101.6, 44.45, 6.35, 76.2, 12.7, 76.2, 76.2, 127, 76.2, 3.18, 6.35, 3.18, 53.98, 9.53, 101.6, 83.90, 50.8, 12.7,
            107.95, 104.14, 35.56, 10.16, 101.6, 38.1, 6.35, 5.08, NA, NA, 7.62, 5.08, 3.18, 5.08, 3.18, 7.62, NA, 12.7, 76.2, 92.20,
            58.67, 44.45, 91.69, 69.85, 6.35, 25.4, 12.7, 69.85, 76.2, 71.12, 6.35, 25.4, 6.35, 71.12, 50.8, 69.85, NA, 27.94, 5.08, 69.85,
            76.2, 60.45, 3.81, 28.7, 7.62, 69.85, 60.33, 69.85, 3.18, 25.4, 6.35, 69.85, 76.2, 76.2, 3.18, 25.4, 6.35, 76.2, 76.2)
actual_a <- c(0.76, 80.01, 3.18, 101.6, 35.56, 6.35, 1.52, 0.76, 0.76, 6.35, 0.76, 1.52, 3.18, 0.76, 6.35, 0.76, 76.2, 80.01, 71.12, 57.15,
              107.95, 0.76, 0.76, 0.76, 0.76, 25.4, 0.76, 80.01, 3.18, 101.6, 35.56, 6.35, 76.2, 107.95, 0.76, 0.76, 0.76, 0.76, 25.4, 76.2,
              57.15, 107.95, 101.6, 35.56, 6.35, 76.2, 12.7, 76.2, 76.2, 127, 76.2, 0.76, 3.18, 0.76, 80.01, 3.18, 76.2, 80.01, 71.12, 57.15,
              107.95, 101.6, 35.56, 6.35, 101.6, 35.56, 6.35, 1.52, 0.76, 0.76, 6.35, 0.76, 1.52, 3.18, 0.76, 6.35, 0.76, 6.35, 76.2, 80.01,
              71.12, 57.15, 107.95, 71.12, 0.76, 35.56, 3.18, 71.12, 76.2, 71.12, 0.76, 35.56, 3.18, 71.12, 76.2, 71.12, 0.76, 35.56, 3.18, 71.12,
              76.2, 71.12, 0.76, 35.56, 3.18, 71.12, 76.2, 71.12, 0.76, 35.56, 3.18, 71.12, 76.2, 71.12, 0.76, 35.56, 3.18, 71.12, 76.2)
censored_a <- c(1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
                1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

data_POD_P12 <- cbind(meas_a,actual_a,censored_a)

POD.calc(data_POD_P12,dist = "Logistic",alth = 0.254,xlabel1 = "Crack Size (mm)",ylabel1 = "POD",unitlabel = "mm")
POD.calc(data_POD_P12,dist = "Loglogistic",alth = 0.254,xlabel1 = "Crack Size (mm)",ylabel1 = "POD",unitlabel = "mm")
POD.calc(data_POD_P12,dist = "Lognormal",alth = 0.254,xlabel1 = "Crack Size (mm)",ylabel1 = "POD",unitlabel = "mm")

# Bayesian Updating
# ===============================================================
library(pracma)
library(StanHeaders)
library(rstan)
library(ggplot2)
library(shinystan)
library(cmdstanr)
library(bayesplot)

# ===========================================================================================
# Log-logistic Inference
# ===========================================================================================
priors <- paste(c("alpha_1 ~ uniform(-10,10);","alpha_2 ~ uniform(-10,10);"),collapse = " ")
POD <- "1./(1 + exp(-(alpha_1 + log(Y).*alpha_2)))"
params <- "real alpha_1; real alpha_2;"
loglik <- paste(c("target += bernoulli_lpmf(Hit|",POD,");"),collapse = "")
paramsvec <- c("alpha_1","alpha_2")

block1 <- "data {int<lower=0> n; array[n] int Hit;  vector[n] Y;}"
datablock <- list(n = length(actual_a), Hit = censored_a, Y = actual_a)
block2 <- paste(c("parameters {",params,"}"),collapse = " ")
block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

stanlscode <- paste(c(block1,block2,block3),collapse=" ")
# stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")

stanlsfile <- write_stan_file(stanlscode)
print(stanlsfile)
# Generate initial list (one list per chain)
pt_est <- c(1,1)
names(pt_est) <- paramsvec
pt_estlist <- as.list(pt_est)
init_pt_est <- vector("list",4)
for(i in 1:4){
  init_pt_est[[i]] <- pt_estlist
}
lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
fit <- sampling(lsmod, data = datablock, iter = 90000, warmup = 1000, init = pt_est)

stats.mean.sd <- summary(fit)$summary[,c(1,3)]
stats.Rhat <- rhat(fit)
confidbounds <- mcmc_intervals_data(data.frame(extract(fit, paramsvec)),prob_outer = 0.9)
outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec),1],unname(stats.mean.sd)[1:length(paramsvec),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec)]), nrow = length(paramsvec), ncol = 6, byrow = FALSE,dimnames = list(paramsvec,c("Mean","Standard Deviation","5%","Median","95%","R\U005E")))
sims <- as.data.frame(fit) # Fit output of Rstan to a data frame for easy modeling

xplot <- linspace(a_lth,max(actual_a),1000)
# predictions <- sapply(1:length(sims$alpha_1), function(i) {
#   1/(1 + exp(-(sims$alpha_1[i] + log(xplot)*sims$alpha_2[i])))
# })
predictions <- sapply(1:10000, function(i) {
  1/(1 + exp(-(sims$alpha_1[i] + log(xplot)*sims$alpha_2[i])))
})
POD.mean <- 1/(1 + exp(-(0.427683 + log(xplot)*(1.612457))))


mean_predictions <- apply(predictions, 1, mean)
lower_ci <- apply(predictions, 1, function(x) quantile(x, 0.05))
upper_ci <- apply(predictions, 1, function(x) quantile(x, 0.95))

df_HitMiss <- data.frame(x = c(actual_a,c(0.76,1.52,3.18,6.35)), y = c(censored_a, c(0.5,0.75,0.92,0.9)),
                         H.or.M = c(rep("Hit or Miss Data",length(actual_a)),rep("POD Point Estimates",4)))
df_POD_curve <- data.frame(x = xplot, y = POD.mean,POD.curve = rep("POD",1000))
df_confidence <- data.frame(x = xplot,ymin = lower_ci, ymax = upper_ci)

plotout <- ggplot() + geom_point(data = df_HitMiss, aes(x, y, shape = H.or.M), size = 2.5, colour = "black") +
  geom_path(data = df_POD_curve, aes(x, y, linetype = POD.curve), linewidth = 0.4, colour = "black") +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),size = .4)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab("Crack Size (mm)") +
  ylab("POD")
plotout <- plotout + geom_ribbon(data = df_confidence, aes(x=x, ymin = ymin, ymax = ymax),alpha = 0.25, fill = "red")

# ===========================================================================================
# Lognormal Inference
# ===========================================================================================
priors <- paste(c("gamma_1 ~ uniform(-10,10);","gamma_2 ~ uniform(0.1,10);"),collapse = " ")
POD <- "0.5.*(1 + erf((log(Y) - gamma_1)./(1.414214.*gamma_2)))"
params <- "real gamma_1; real gamma_2;"
loglik <- paste(c("target += bernoulli_lpmf(Hit|",POD,");"),collapse = "")
paramsvec <- c("gamma_1","gamma_2")

block1 <- "data {int<lower=0> n; array[n] int Hit;  vector[n] Y;}"
datablock <- list(n = length(actual_a), Hit = censored_a, Y = actual_a)
block2 <- paste(c("parameters {",params,"}"),collapse = " ")
block3 <- paste(c("model {",priors,loglik,"}"),collapse = " ")

stanlscode <- paste(c(block1,block2,block3),collapse=" ")
# stanlscode <- paste(c(block1,block2,block2b,block3),collapse=" ")

stanlsfile <- write_stan_file(stanlscode)
print(stanlsfile)
# Generate initial list (one list per chain)
pt_est <- c(1,1)
names(pt_est) <- paramsvec
pt_estlist <- as.list(pt_est)
init_pt_est <- vector("list",4)
for(i in 1:4){
  init_pt_est[[i]] <- pt_estlist
}
lsmod <- stan_model(model_code = stanlscode, verbose = TRUE)
fit <- sampling(lsmod, data = datablock, iter = 90000, warmup = 1000, init = pt_est)

stats.mean.sd <- summary(fit)$summary[,c(1,3)]
stats.Rhat <- rhat(fit)
confidbounds <- mcmc_intervals_data(data.frame(extract(fit, paramsvec)),prob_outer = 0.9)
outputtable <- matrix(c(unname(stats.mean.sd)[1:length(paramsvec),1],unname(stats.mean.sd)[1:length(paramsvec),2],confidbounds[[5]],confidbounds[[7]],confidbounds[[9]],unname(stats.Rhat)[1:length(paramsvec)]), nrow = length(paramsvec), ncol = 6, byrow = FALSE,dimnames = list(paramsvec,c("Mean","Standard Deviation","5%","Median","95%","R\U005E")))
sims <- as.data.frame(fit) # Fit output of Rstan to a data frame for easy modeling

xplot <- linspace(a_lth,max(actual_a),1000)
# predictions <- sapply(1:length(sims$alpha_1), function(i) {
#   1/(1 + exp(-(sims$alpha_1[i] + log(xplot)*sims$alpha_2[i])))
# })
predictions <- sapply(1:10000, function(i) {
  plnorm(xplot,sims$gamma_1[i],sims$gamma_2[i])
})
POD.mean <- plnorm(xplot,-0.5262002,1.6433772)


mean_predictions <- apply(predictions, 1, mean)
lower_ci <- apply(predictions, 1, function(x) quantile(x, 0.05))
upper_ci <- apply(predictions, 1, function(x) quantile(x, 0.95))

df_HitMiss <- data.frame(x = c(actual_a,c(0.76,1.52,3.18,6.35)), y = c(censored_a, c(0.5,0.75,0.92,0.9)),
                         H.or.M = c(rep("Hit or Miss Data",length(actual_a)),rep("POD Point Estimates",4)))
df_POD_curve <- data.frame(x = xplot, y = POD.mean,POD.curve = rep("POD",1000))
df_confidence <- data.frame(x = xplot,ymin = lower_ci, ymax = upper_ci)

plotout2 <- ggplot() + geom_point(data = df_HitMiss, aes(x, y, shape = H.or.M), size = 2.5, colour = "black") +
  geom_path(data = df_POD_curve, aes(x, y, linetype = POD.curve), linewidth = 0.4, colour = "black") +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),size = .4)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab("Crack Size (mm)") +
  ylab("POD")
plotout2 <- plotout2 + geom_ribbon(data = df_confidence, aes(x=x, ymin = ymin, ymax = ymax),alpha = 0.25, fill = "red")
