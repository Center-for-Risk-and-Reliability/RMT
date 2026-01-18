# Distribution Rank and Fit
# Developed by Dr. Reuel Smith, 2025

distribution.fit <- function(x){
  data <- cbind(x,rep(1,length(x)),rep(1,length(x))) # Set up data table for analysis
  OUT1 <- probplot.wbl(data,"Blom",xlabel1 = "x",MLE_i = 1)   # 1. Weibull
  OUT2 <- probplot.wbl3P(data,"Blom",xlabel1 = "x",MLE_i = 1) # 2. Three Parameter Weibull
  OUT3 <- probplot.exp(data,"Blom",xlabel1 = "x",MLE_i = 1) # 3. Exponential
  OUT4 <- probplot.exp2P(data,"Blom",xlabel1 = "x",MLE_i = 1) # 4. Two Parameter Exponential
  OUT5 <- probplot.nor(data,"Blom",xlabel1 = "x",MLE_i = 1) # 5. Normal
  OUT6 <- probplot.logn(data,"Blom",xlabel1 = "x",MLE_i = 1) # 6. Lognormal
  OUT7 <- probplot.gam(data,"Blom",xlabel1 = "x",MLE_i = 1) # 7. Gamma
  OUT8 <- probplot.gam3P(data,"Blom",xlabel1 = "x",MLE_i = 1) # 8. Three Parameter or Generalized Gamma
  OUT9 <- probplot.logist(data,"Blom",xlabel1 = "x",MLE_i = 1) # 9. Logistic
  OUT10 <- probplot.loglogist(data,"Blom",xlabel1 = "x",MLE_i = 1) # 10. Log-logistic
  OUT11 <- probplot.gumb(data,"Blom",xlabel1 = "x",MLE_i = 1) # 11. Gumbel

  parameters <- list(Weibull.MLE = OUT1$output[[1]]$`Parameter Estimates`,
                     Three.Parameter.Weibull.MLE = OUT2$output[[1]]$`Parameter Estimates`,
                     Exponential.MLE = OUT3$output[[1]]$`Parameter Estimates`,
                     Two.Parameter.Exponential.MLE = OUT4$output[[1]]$`Parameter Estimates`,
                     Normal.MLE = OUT5$output[[1]]$`Parameter Estimates`,
                     Lognormal.MLE = OUT6$output[[1]]$`Parameter Estimates`,
                     Gamma.MLE = OUT7$output[[1]]$`Parameter Estimates`,
                     Generalized.Gamma.MLE = OUT8$output[[1]]$`Parameter Estimates`,
                     Logistic.MLE = OUT9$output[[1]]$`Parameter Estimates`,
                     Loglogistic.MLE = OUT10$output[[1]]$`Parameter Estimates`,
                     Gumbel.MLE = OUT11$output[[1]]$`Parameter Estimates`)
  loglikelihood <- list(Weibull.loglik = OUT1$output[[1]]$loglikelihood,
                        Three.Parameter.Weibull.loglik = OUT2$output[[1]]$loglikelihood,
                        Exponential.loglik = OUT3$output[[1]]$loglikelihood,
                        Two.Parameter.Exponential.loglik = OUT4$output[[1]]$loglikelihood,
                        Normal.loglik = OUT5$output[[1]]$loglikelihood,
                        Lognormal.loglik = OUT6$output[[1]]$loglikelihood,
                        Gamma.loglik = OUT7$output[[1]]$loglikelihood,
                        Generalized.Gamma.loglik = OUT8$output[[1]]$loglikelihood,
                        Logistic.loglik = OUT9$output[[1]]$loglikelihood,
                        Loglogistic.loglik = OUT10$output[[1]]$loglikelihood,
                        Gumbel.loglik = OUT11$output[[1]]$loglikelihood)
  AIC <- list(Weibull.loglik = OUT1$output[[1]]$AIC,
                        Three.Parameter.Weibull.loglik = OUT2$output[[1]]$AIC,
                        Exponential.loglik = OUT3$output[[1]]$AIC,
                        Two.Parameter.Exponential.loglik = OUT4$output[[1]]$AIC,
                        Normal.loglik = OUT5$output[[1]]$AIC,
                        Lognormal.loglik = OUT6$output[[1]]$AIC,
                        Gamma.loglik = OUT7$output[[1]]$AIC,
                        Generalized.Gamma.loglik = OUT8$output[[1]]$AIC,
                        Logistic.loglik = OUT9$output[[1]]$AIC,
                        Loglogistic.loglik = OUT10$output[[1]]$AIC,
                        Gumbel.loglik = OUT11$output[[1]]$AIC)

  return(list(parameters=parameters,loglikelihood,AIC=AIC))

}
