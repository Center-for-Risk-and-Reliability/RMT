# Least-Squares PH-Weibull Life-Stress Estimator
# Developed by Dr. Reuel Smith, 2022

lifestress.PHWbl.LSQest <- function(pp,interact_stress){
  #Load pracma library for pseudoinverse
  library(pracma)

  # First check and see that there are multiple stress levels
  if(length(pp)<3) {
    stop('Need more than one stress level to generate estimates')
  }
  # Then check and see if there are single entry data
  if(length(pp)%%3==0){
    singledat<-0 # FALSE Single data does not exist
  } else{
    singledat<-1 # TRUE Single data exists
  }

  # Sets the default computation for stresses to simply single evaluation
  # and not check for interactions
  if(missing(interact_stress)){
    interact_stress<-0
  }

  # Setup vectors (for cases with and without single point data)
  if(singledat==0){
    # Sets up existing probability plot curve life and stress vectors
    L<-rep(0,length(pp)/3)
    if (length(pp[[1]])<2){
      S<-rep(0,length(pp)/3)
    } else {
      S<-matrix(rep(0,(length(pp)/3)*length(pp[[1]])),nrow=length(pp)/3,ncol=length(pp[[1]]),byrow = TRUE)
    }
    distparams<-rep(0,length(pp)/3)
  } else if(singledat==1){
    # Sets up probability plot curve and single entry L-S life and stress vectors
    L<-rep(0,(length(pp)-1)/3 + length(tail(pp,n=1)[[1]]))
    if (length(pp[[1]])<2){
      S<-rep(0,(length(pp)-1)/3 + length(tail(pp,n=1)[[1]]))
    } else {
      # NOTE TEST THIS UNDER APPROPRIATE CIRCUMSTANCES
      S<-matrix(rep(0,((length(pp)-1)/3 + length(tail(pp,n=1)[[1]]))*length(pp[[1]])),nrow=(length(pp)-1)/3 + length(tail(pp,n=1)[[1]]),ncol=length(pp[[1]]),byrow = TRUE)
    }
    # Distribution parameter pulls ONLY apply to the probability plots
    distparams<-rep(0,(length(pp)-1)/3)
  }


  # Fill in Stress and Life Vectors
  if(singledat==0){
    for(i2 in 1:(length(pp)/3)){
      # Stress Levels
      if (length(pp[[1]])<2){
        S[i2]<-pp[[i2*3-2]]
      } else {
        for(j in 1:length(pp[[1]])){
          S[i2,j] <- pp[[i2*3-2]][[j]]
        }
      }

      # Life Estimates
      L[i2]<-pp[[i2*3-1]][,1]
      distparams[i2]<-pp[[i2*3-1]][2]
    }
  } else if(singledat==1){
    # First Tabulate Probability Plot S and L data
    for(i2 in 1:((length(pp)-1)/3)){
      # Stress Levels
      if (length(pp[[1]])<2){
        S[i2]<-pp[[i2*3-2]]
      } else {
        for(j in 1:length(pp[[1]])){
          S[i2,j] <- pp[[i2*3-2]][[j]]
        }
      }
      # Life Estimates
      L[i2]<-pp[[i2*3-1]][,1]
      distparams[i2]<-pp[[i2*3-1]][2]
    }
    # Next tabulate the single point data
    for(i2 in 1:length(tail(pp,n=1)[[1]])){
      S[i2+(length(pp)-1)/3]<-tail(pp,n=1)[[1]][[i2]][,3]
      L[i2+(length(pp)-1)/3]<-tail(pp,n=1)[[1]][[i2]][,1]
    }
  }

  # Writeup for the output text
  ls_txt<-"Parametric Proportional Hazard"
  pdf_txt<-"\U03B2\U2219 t^(\U03B2-1)\U2219 t\U2219 exp[\U03A3_(j=0)^m \U03B8_j S_j - t^\U03B2 exp(\U03A3_(j=0)^m \U03B8_j S_j)]"
  dist_txt<-"Weibull"
  distparam_txt<-"\U03B2"

  # Set up Matrix and Vector to perform pseudo-inverse operation to get LSQ parameters
  # theta_0, ... theta_n
  # Main vector -beta*LN(alpha)
  V_main <- -distparams*log(L)

  stress_txt<-paste("\U03B8", 0:length(pp[[1]]), sep = "_")

  if(interact_stress==0){
    # No Stress dependency check
    M_main <- matrix(c(rep(1,length(pp)/3),c(S)), nrow = length(pp)/3, ncol = 1+length(pp[[1]]), byrow = FALSE)
  }
  if(interact_stress==1){
    # Stress dependency check
    # Permutations of stress pairs
    Spairs <- combn(length(pp[[1]]),2)
    for(i in 1:dim(Spairs)[2]){
      if(i==1){
        Spair_v <- S[, Spairs[1,i]]*S[, Spairs[2,i]]
        stressjoint_txt<-paste("\U03B8", paste(Spairs[,i],collapse=""), sep = "_")
      } else{
        Spair_v <- c(Spair_v,S[, Spairs[1,i]]*S[, Spairs[2,i]])
        stressjoint_txt<-c(stressjoint_txt,paste("\U03B8", paste(Spairs[,i],collapse=""), sep = "_"))
      }
    }
    M_main <- matrix(c(rep(1,length(pp)/3),c(S),Spair_v), nrow = length(pp)/3, ncol = 1+length(pp[[1]])+dim(Spairs)[2], byrow = FALSE)
    stress_txt<-c(stress_txt,stressjoint_txt)
    S<-M_main[,2:(1+length(pp[[1]])+dim(Spairs)[2])]
  }

  # Compute the theta LSQ parameters
  params  <- pinv(M_main)%*%V_main
  LSQ <- c(mean(distparams),params)

  params_txt<-c(distparam_txt,stress_txt)

  # Produce some output text that summariZes the results
  cat(c("Least-Squares estimates for the ",ls_txt,"-",dist_txt," Life-Stress model.\n\nf(t|S;\U03B8) = ",pdf_txt,"\n\n"),sep = "")
  print(matrix(c(LSQ), nrow = 1, ncol = length(LSQ), byrow = TRUE,dimnames = list(c("Life-Stress Parameters"),params_txt)))
  cat("\n\n")

  # Return parameter list
  return(list(S,L,LSQ))
  # return(list(S,L,distparams,V_main,M_main))
}
