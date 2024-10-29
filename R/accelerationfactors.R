# Acceleration Factor Calculator
# Developed by Dr. Reuel Smith, 2021-2024

accelfactor <- function(params,lsm,Sacc,Suse){
  lsfcn<-lifestress.select(lsm)[[1]]
  if(length(Suse)==1 || lsm=="MultiStress") {
    if(length(Sacc) == 1){
      AF<-lsfcn(params,Suse)/lsfcn(params,Sacc)
    } else{
      AF <- rep(0,length(Sacc))
      for(i in 1:length(Sacc)){
        AF[i] <- lsfcn(params,Suse)/lsfcn(params,Sacc[i])
      }
    }
  }
  if(length(Suse)==2 && !(lsm=="MultiStress")) {
    if(is.null(dim(Sacc)) == TRUE){
      AF<-lsfcn(params,Suse[1],Suse[2])/lsfcn(params,Sacc[1],Sacc[2])
    }
    if(is.null(dim(Sacc)) == FALSE){
      AF<-rep(0,length(Sacc[,1]))
      for(i in 1:length(Sacc[,1])){
        AF[i] <- lsfcn(params,Suse[1],Suse[2])/lsfcn(params,Sacc[i,1],Sacc[i,2])
      }
    }
  }
  return(AF)
}
