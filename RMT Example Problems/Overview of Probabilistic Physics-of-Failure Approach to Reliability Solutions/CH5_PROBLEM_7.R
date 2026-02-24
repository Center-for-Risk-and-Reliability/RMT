# CHAPTER 5 PROBLEM 7
# Reuel Smith
# =================================================
# PART a)
# =================================================================
# RMT Code Entry
# =================================================================
data_ADT_PROBLEM_7 <-cbind(data.frame(time=rep(c(10, 20, 30, 40, 50, 60, 70),5),
                                      deg=c(240.7, 316.2, 380.7, 409.1, 419.2, 455.8, 475.8,
                                            255.5, 297.9, 392.4, 436, 459.8, 472.8, 479.8,
                                            293.3, 316.1, 343.7, 397.5, 429.2, 444, 449.3,
                                            268.7, 305, 338.2, 359.5, 442.5, 464.6, 468.9,
                                            295.5, 348.1, 362, 400.9, 445.3, 455.3, 470.2),
                                      unitno=c(rep('Switch 1',7),rep('Switch 2',7),rep('Switch 3',7),rep('Switch 4',7),rep('Switch 5',7)),
                                      stress=c(rep(1,35))))

adt.rank(data_ADT_PROBLEM_7)   # Run the ADT ranking tool to determine which of the now eight degradation-life models is the better fit
Df.PROBLEM_7 <- 550*0.9        # Endurance limit or resistance in ohms (495 ohms)
# Least Squares evaluation
degradationlife.LSQest(data=data_ADT_PROBLEM_7,dl="SquareRoot2",dist="Normal",
                       pp="Blom",Df=495, modelstress = NULL,
                       xlabel = "Time (hours)",ylabel = "Reistance (ohms)")
# MLE evaluation
degradationlife.MLEest(data=data_ADT_PROBLEM_7,dl="SquareRoot2",dist="Normal",
                       pp="Blom",Df=495, modelstress = NULL,confid=0.95,
                       xlabel = "Time (hours)",ylabel = "Reistance (ohms)")

# PART b
# =================================================================

# Simulation by MLE
library(mvtnorm)                                # Load mvtnorm to run rmvnorm function
Output.PROBLEM_7 <- degradationlife.MLEest(data=data_ADT_PROBLEM_7,dl="SquareRoot2",dist="Normal",
                                  pp="Blom",Df=495, modelstress = NULL,confid=0.95,
                                  xlabel = "Time (hours)",ylabel = "Reistance (ohms)")
theta.hat <- Output.PROBLEM_7$MLE.point.estimate
Var.cov <- Output.PROBLEM_7$var.cov.matrix      # Extract the variance/covariance matrix
N <- 1000                                       # State sample number that will be extracted
param_pull <- rmvnorm(N,theta.hat,Var.cov)
fail_times <- rep(0,N)

for(i in 1:N){
  # Compute failure times based on each parameter pull
  fail_times[i] <- ((495 - param_pull[i,2])/param_pull[i,3])^2
}
# Check to see best distribution to fit time-to-failure
probplot.nor(cbind(fail_times,rep(1,1000),rep(1,1000)),pp="Blom",xlabel1 = "Time to Failure (hours)",MLE_i = 1)
# $output[[1]]$loglikelihood
# [1] -2307.709
probplot.logn(cbind(fail_times,rep(1,1000),rep(1,1000)),pp="Blom",xlabel1 = "Time to Failure (hours)",MLE_i = 1)
# $output[[1]]$loglikelihood
# [1] -2305.158     <======== HIGHEST OF THREE
probplot.wbl(cbind(fail_times,rep(1,1000),rep(1,1000)),pp="Blom",xlabel1 = "Time to Failure (hours)",MLE_i = 1)
# $output[[1]]$loglikelihood
# [1] -2396.827

plnorm(155, 4.335971, 0.03175352)
