# Short Script for Chapter 5 Problem 6
# Reuel Smith

library(nls.multstart)
P2MLE <- adt.full.MLE(data641FinalP2,"CrackProp2","Lognormal",0.034663947,0.95,"twosided")
Nf_samples2 <- linspace(100000,800000,10000)
R_deg2 <- 1 - plnorm(0.034663947,log(0.001) + (pi*exp(P2MLE[[1]][1])*(200^2)*Nf_samples2),exp(P2MLE[[1]][2]))
F_Nf <- R_deg2

Nf_dist <- function(mut,sigt){
  plnorm(Nf_samples2,mut,sigt)
}
df2 <- data.frame(y = F_Nf)
dist_est <- nls_multstart(y ~ Nf_dist(mut,sigt), data = df2, iter = 5000, start_lower = c(mut = 9,sigt = 0),start_upper = c(mut = 20,sigt = 10), lower = c(mut = -100,sigt = 0))
dist_est$m$getPars()
