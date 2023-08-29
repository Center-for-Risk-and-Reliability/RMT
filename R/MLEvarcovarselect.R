# Maximum Likelihood Variance-Covariance Matrix Select
# Developed by Dr. Reuel Smith, 2022

MLE.var.covar.select <- function(loglik,LSQest){
  library(pracma)
  library(matrixcalc)
  library(ucminf)

  # There is still some question as to which R script and library to use for MLE estimation.
  # nlm provides the mean estimate at all times, however there are cases that produce a negative
  # variance in some variables.  Further study will need to be done to find the best option.  We may
  # need to have multiple options on standby just in case.  There is the 'hessian' function from
  # the numDeriv library as a first alternative.  A second alternative is the 'ucminf' function
  # from the library of the same name
  # fish <- hessian(loglik,out$estimate)

  # TIER 1 ESTIMATE: Based on nlm
  out1 <- nlm(loglik, theta <- LSQest, hessian=TRUE)
  theta.hat1 <- out1$estimate
  fish1 <- out1$hessian
  inv.fish1 <- pinv(fish1)

  if(min(diag(inv.fish1)) > 0){
    inv.fish <- inv.fish1
    theta.hat <- theta.hat1}
  else {
    # TIER 2 ESTIMATE: Based on ucminf
    out2 <- ucminf(LSQest,loglik,hessian=1)
    theta.hat2 <- out2$par
    fish2 <- out2$hessian
    inv.fish2 <- pinv(fish2)

    if(min(diag(inv.fish2)) > 0){
      inv.fish <- inv.fish2
      theta.hat <- theta.hat2
    } else{
      # TIER 3 ESTIMATE: Based on optim
      out3 <- optim(LSQest,loglik,hessian=1)
      theta.hat3 <- out3$par
      fish3 <- out3$hessian
      inv.fish3 <- pinv(fish3)

      if(min(diag(inv.fish3)) > 0){
        inv.fish <- inv.fish3
        theta.hat <- theta.hat3
      } else {
        inv.fish <- inv.fish1
        theta.hat <- theta.hat1
      }
    }
  }

  # ESTIMATE 2: Based on ucminf
  # out2 <- ucminf(LSQest,loglik,hessian=1)
  # theta.hat2 <- out2$par
  # fish2 <- out2$hessian
  # inv.fish2 <- pinv(fish2)
  #
  # # ESTIMATE 3: Based on optim
  # out3 <- optim(LSQest,loglik,hessian=1)
  # theta.hat3 <- out3$par
  # fish3 <- out3$hessian
  # inv.fish3 <- pinv(fish3)
  #
  # # Pick the best estimate to work with (will add more options as I find them)
  # if(min(diag(inv.fish1)) > 0){
  #   inv.fish <- inv.fish1
  #   theta.hat <- theta.hat1
  # } else if(min(diag(inv.fish2)) > 0){
  #   inv.fish <- inv.fish2
  #   theta.hat <- theta.hat2
  # } else if(min(diag(inv.fish3)) > 0){
  #   inv.fish <- inv.fish3
  #   theta.hat <- theta.hat3
  # } else {
  #   inv.fish <- inv.fish1
  #   theta.hat <- theta.hat1
  # }

  return(list(theta.hat,inv.fish))
}
