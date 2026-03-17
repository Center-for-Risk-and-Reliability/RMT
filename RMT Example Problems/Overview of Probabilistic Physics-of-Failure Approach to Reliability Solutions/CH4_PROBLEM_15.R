# CHAPTER 4 PROBLEM 15
# Gabriel Cruz
# =================================================
# MLE Part a #
# Initialize data and constants #
# Definitions of different variables used in this code
# x[1]: loga parameter
# x[2]: n parameter
# x[3]: sigma_t parameter
# Stress in MPa V <- c(375*(c(1:7)/c(1:7)),425*(c(1:7)/c(1:7)))
# The given ungrouped failure times
ttf <- c(32,36,40,43,48,54,58,16,18,18,20,22,26,28)
mtot <- length(V)
loglik <- function(X) {-(-mtot*log(X[3]) - mtot*0.5*log(2*pi) - sum(log(ttf)) - 0.5*sum(((log(ttf) + log(exp(X[1])) + X[2]*log(V))/X[3])^2)) }
output <- nlm(loglik, theta <- c(-20,3,0.5), hessian=TRUE)
est <- output$estimate
fish <- output$hessian
a_1 <- exp(est[1])
n_1 <- est[2]
sigmat_1 <- est[3]
Nat300_1 <- exp(-log(a_1) - n_1*log(300))
Part b).
# SETUP FOR SAMPLING #
niter <- 600000
# Parameter 1 Priors: a_L and a_H. Parameter 1 A is assumed to have a # uniform prior distribution UNIF(loga_L,loga_H).
#logA_L <- -20
logA_L <- log(1e-19)
logA_H <- -16
# Parameter 2 Priors: B_mu and B_tau. Parameter 2 B is assumed to have a
# normal prior distribution
NORM(B_mu,B_tau) where B_tau = 1/B_SD^2.
#n_L <- 3
n_L <- 0.3
n_H <- 8
# Parameter 3 Priors: beta_L and beta_H. Parameter 3 beta is assumed to have a
# normal prior distribution
NORM(beta_mu,beta_tau) where beta_tau = 1/beta_SD^2.
sigma_L <- 0.1
sigma_H <- 0.4
t <- ttf
nt <- length(t)
S <- V
data <- list("logA_L","logA_H","n_L","n_H","sigma_L","sigma_H","t","nt","S")
parameters <- c("logA","n","sigma")
inits <- list(list(logA = -18, n = 5, sigma = 0.2))
# Under MODEL.FILE, change the folder name to where your WinBUGS model
# is located. Under BUG.DIRECTORY, change the folder name to where the
# WinBUGS executable is located.
ADT.sim <- bugs(data,inits,parameters,model.file="C:/…",n.chains = 1, n.iter = niter, n.burnin=25000, n.thin=16, debug=TRUE, bugs.directory="C:/… ", DIC=FALSE)
attach.bugs(ADT.sim)
#detach.bugs()
flush.console()
a_2 <- exp(logA)
n_2 <- n
sigmat_2 <- sigma
Nat300_2 <- exp(-log(a_2) - n_2*log(300))
conf <- 0.9
lowlim <- (1-conf)/2
highlim <- 1-lowlim
results1 <- matrix(c(quantile(a_2,lowlim),median(a_2),quantile(a_2,highlim),mean(a_2),sd(a_2),quantile(n_2,lowlim ),median(n_2),quantile(n_2,highlim),mean(n_2),sd(n_2),quantile(sigmat_2,lowlim),median(sigmat_2),qu antile(sigmat_2,highlim),mean(sigmat_2),sd(sigmat_2),quantile(Nat300_2,lowlim),median(Nat300_2),qu antile(Nat300_2,highlim),mean(Nat300_2),sd(Nat300_2)), nrow = 4, ncol = 5, byrow = TRUE,dimnames = list(c("a", "n", "sigma_t", "Life at 300 MPa"),c("5%", "median", "95%", "mean", "SD")))
