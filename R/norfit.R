# Normal Fit
# Developed by Dr. Reuel Smith, 2023

fit.nor <- function(xi, rc = NULL, pp = "Blom", confid = 0.95, sided = "twosided"){
  library(pracma)

  conf.level <- confid

  # Initial LSQ estimate
  x_R <- plotposit.select(xi, rc, pp)
  est1 <- c(probplotparam.nor(x_R[,1], x_R[,2])[[3]])
  # MLE estimate
  est2 <- distribution.MLEest(est1, "Normal", xi, rc, confid, sided)
  # AIC and BIC estimates
  AIC = 2*2 - 2*(est2[[4]])
  BIC = 2*log(length(xi)+length(rc)) - 2*(est2[[4]])

  if(sided == "twosided"){
    CIlow <- c(est2[[3]][[1]][1],est2[[3]][[2]][1])
    CIhigh <- c(est2[[3]][[1]][2],est2[[3]][[2]][2])
    conflim_txt<-c(paste(c("Lower ",100*conf.level,"%"),collapse = ""),paste(c("Upper ",100*conf.level,"%"),collapse = ""))
  }
  if(sided == "onesidedlow"){
    CIone <- c(est2[[3]][[1]],est2[[3]][[2]])
    conflim_txt<-paste(c("One-Sided Low ",100*conf.level,"%"),collapse = "")
  }
  if(sided == "onesidedhigh"){
    CIone <- c(est2[[3]][[1]],est2[[3]][[2]])
    conflim_txt<-paste(c("One-Sided High ",100*conf.level,"%"),collapse = "")
  }
  # Produce some output text that summarizes the results
  cat(c("Maximum-Likelihood estimates for Normal Fit.\n\n"),sep = "")
  if(sided == "twosided"){
    print(matrix(c(est2[[1]],sqrt(diag(est2[[2]])/(length(xi)+length(rc))),CIlow,CIhigh), nrow = 4, ncol = 2, byrow = TRUE,dimnames = list(c("Parameters Point Estimate","Standard Error",conflim_txt),c("\U03BC (Location)","\U03C3 (Scale)"))))
  } else{
    print(matrix(c(est2[[1]],sqrt(diag(est2[[2]])/(length(xi)+length(rc))),CIone), nrow = 3, ncol = 2, byrow = TRUE,dimnames = list(c("Parameters Point Estimate","Standard Error",conflim_txt),c("\U03BC (Location)","\U03C3 (Scale)"))))
  }
  cat(c("\n"),sep = "")


  return(list(normparams = est2[[1]], CImu = est2[[3]][[1]], CIsigma = est2[[3]][[2]], varcovmat = est2[[2]], loglikelihood = est2[[4]], AIC = AIC, BIC = BIC))
}
