# CHAPTER 5 PROBLEM 5
# Coded by Dr. Reuel Smith 2019
# =================================================

# ln(β) - parameter 1 theta[1]
# α - parameter 2 theta[2]
# sig - parameter 3 theta[3]

# Data Entry
TTF <- c(10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000)
r <- c(0.01, 0.02, 0.05, 0.07, 0.09, 0.1, 0.005, 0.017, 0.018, 0.03, 0.04, 0.065, 0.02, 0.04, 0.055, 0.065, 0.078, 0.09, 0.018, 0.035, 0.04, 0.062, 0.087, 0.195, 0.03, 0.05, 0.065, 0.089, 0.12, 0.15)
Temp <- c(395, 395, 395, 395, 395, 395, 390, 390, 390, 390, 390, 390, 500, 500, 500, 500, 500, 500, 385, 385, 385, 385, 385, 385, 390, 390, 390, 390, 390, 390)
pH <- c(5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5)

loglik <- function(theta) {
  -sum(-log(r) - log(theta[3]) - 0.5*log(2*pi) - 0.5*(1/theta[3]^2)*(log(r) - (1/3)*(log(TTF) - theta[1] - log(pH) -theta[2]/Temp))^2)
}

guess <- c(13.555,0.031,1)
nlm(loglik, theta <- guess, hessian=TRUE)
out <- nlm(loglik, theta <- guess, hessian=TRUE)
fish <- out$hessian
solve(fish)

theta.hat <- out$estimate
theta.hat
conf.level <- 0.95
crit <- qnorm((1 + conf.level)/2)
inv.fish <- solve(fish)

# Upper and lower limits of the parameters

theta.hat[1] + c(-1, 1) * crit * sqrt(inv.fish[1, 1])
theta.hat[2] + c(-1, 1) * crit * sqrt(inv.fish[2, 2])
theta.hat[3] + c(-1, 1) * crit * sqrt(inv.fish[3, 3])

# =================================================================
# RMT Code Entry
# =================================================================
data_ADT_PROBLEM_5 <-cbind(data.frame(time=rep(c(10, 100, 150, 200, 500, 1000),5),
                                deg=c(0.01, 0.02, 0.05, 0.07, 0.09, 0.1,
                                      0.005, 0.017, 0.018, 0.03, 0.04, 0.065,
                                      0.02, 0.04, 0.055, 0.065, 0.078, 0.09,
                                      0.018, 0.035, 0.04, 0.062, 0.087, 0.195,
                                      0.03, 0.05, 0.065, 0.089, 0.12, 0.15),
                                unitno=c(rep('Pit 1',6),rep('Pit 2',6),rep('Pit 3',6),rep('Pit 4',6),rep('Pit 5',6)),
                                temp=c(rep(395,6),rep(390,6),rep(500,6),rep(385,6),rep(390,6)),
                                PH=c(rep(5.5,6),rep(5,6),rep(4,6),rep(4.5,6),rep(3.5,6))))
degradationlife.LSQest(data=data_ADT_PROBLEM_5,dl="KondoWei",dist="Normal",
                       pp="Blom",Df=0.05, modelstress = NULL,
                       xlabel = "Time (hours)",ylabel = "Pit semi-spherical radius (mm)")

# Short Code for Part c
library(pracma)
N <- 1000
L.pull <- exp(rnorm(N,11.424402,3.11963366))*6.95*exp(rnorm(N,570.2725,1267.160577)/270)
t5 <- quantile(L.pull,0.05)

