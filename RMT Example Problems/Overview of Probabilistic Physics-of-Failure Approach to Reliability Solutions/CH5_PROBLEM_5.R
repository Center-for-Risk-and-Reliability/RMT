# MLE Degradation Problem
# Coded by Dr. Reuel Smith 2019

# ln(A) - parameter 1 theta[1]
# B - parameter 2 theta[2]
# sig - parameter 3 theta[3]

TTF <- c(10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000, 10, 100, 150, 200, 500, 1000)
r <- c(0.01, 0.02, 0.05, 0.07, 0.09, 0.1, 0.005, 0.017, 0.018, 0.03, 0.04, 0.065, 0.02, 0.04, 0.055, 0.065, 0.078, 0.09, 0.018, 0.035, 0.04, 0.062, 0.087, 0.195, 0.03, 0.05, 0.065, 0.089, 0.12, 0.15)
pH <- c(5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5)
T <- c(395, 395, 395, 395, 395, 395, 390, 390, 390, 390, 390, 390, 500, 500, 500, 500, 500, 500, 385, 385, 385, 385, 385, 385, 390, 390, 390, 390, 390, 390)

loglik <- function(theta) {
  -sum(-log(r) - log(theta[3]) - 0.5*log(2*pi) - 0.5*(1/theta[3]^2)*(log(r) - (1/3)*(log(TTF) - theta[1] - log(pH) -theta[2]/T))^2)
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
