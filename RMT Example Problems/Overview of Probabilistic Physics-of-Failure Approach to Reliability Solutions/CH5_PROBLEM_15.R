# CHAPTER 5 PROBLEM 15
# Reuel Smith
# =================================================
library(reliabilityRMT)

# PART 1
a_meas <- c(0.5, 1, 1, 1, 1, 1, 1, 1.1, 1.2, 1.3,
            1.3, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
            1.5, 1.5, 1.5, 1.5, 1.75, 1.8, 1, 2, 2,
            2, 2, 2, 1.1, 2.2, 2.3, 2.5, 3, 3, 3, 3,
            3.5, 4, 3.5, 4.5)

# Update the knowledge of λ hyper-parameter with new data
out_CH5_P15 <- distribution.BAYESest(pt_est = 0.4, dist = "Exponential", TTF = a_meas, Tc = NULL,Tlc=NULL,
                                     confid = 0.95, priors = c("gamma(8.9,990)"),
                                     nsamples = 50000,burnin = 1000,nchains = 4)
Post.dat <- extract(out_CH5_P15$posterior.fit,c("lambda"))$lambda
distribution.fit(Post.dat)

#                     alpha         beta
# Gamma Parameters 51.21455 0.0009311599
# [[2]]$Gamma.loglik
# [1] 705338.8

priorset.2 <- c("gamma(8.9, 0.001010101)")
Dat.out.2<-distribution.BAYESest(pt_est = c(0.05),
                                 "Exponential",TTF= c(430,534,560,560,1403,2020),Tc = rep(2160,4),
                                 confid = 0.9,
                                 priors = priorset.2,
                                 nsamples = 20000,
                                 burnin = 1000,nchains = 4)
posterior_ver.1 <- density(c(Post.dat))
# prior_ver.1 <- dgamma(posterior_ver.1$x,8.9,1/0.001010101)
prior_ver.1 <- density(rgamma(length(posterior_ver.1$x),8.9,1/0.001010101))
df.posterior.prior <- data.frame(x = c(prior_ver.1$x,posterior_ver.1$x),
                                 ymin = c(rep(0,length(prior_ver.1$x)),rep(0,length(posterior_ver.1$x))),
                                 ymax = c(prior_ver.1$y,posterior_ver.1$y),
                                 distlabel = c(rep("Prior",length(prior_ver.1$x)),rep("Posterior",length(posterior_ver.1$x))))
plot.P5.15.Bayes.compare <- ggplot() + geom_ribbon(data = df.posterior.prior, aes(x=x, ymin = ymin, ymax = ymax, fill = distlabel) ,alpha = 0.5) +
  theme(panel.background = element_rect(fill = NA),panel.grid = element_line(colour = "grey80"),axis.line = element_line(arrow = arrow(length = unit(0.05, "inches")),linewidth = .4)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  xlab("λ - hyper parameter") +
  ylab("density")
plot.P5.15.Bayes.compare

# ================================================================================
# PART 2
n_hat <- c(0,19,rep(0,7))                              # observed flaws
a_low <- c(0,1.9,3.2,4.4,5.7,7.0,8.3,9.5,10.8)         # low range of flaw sizes (mm)
a_high <- c(1.9,3.2,4.4,5.7,7.0,8.3,9.5,10.8,12.1)     # upper range of flaw sizes (mm)
a_meas <- a_low + 0.5*(a_high - a_low)                 # measured flaw size (median)
a_true <- rep(0,length(n_hat))
# Estimate true # of flaws
n_true_est <- rep(0,length(n_hat))
for (i in 1:length(n_hat)){
  while(a_true[i] <= a_low[i] || a_true[i] >= a_high[i]){
    a_true[i] <- 1.192*a_meas[i] - 1.08 + rnorm(1,0,0.757)  # Estimate for true flaw size
  }
  if(n_hat[i] > 0){
    n_true_est[i] <- n_hat[i]/(1/(1 + exp(-((a_true[i] - 0.1124)/0.508)))) # estimate for actual number of flaws
  } else {
    n_true_est[i] <- 0
  }
}
n_true_est <- ceil(n_true_est)                          # estimate for actual number of flaws
n_true_est

