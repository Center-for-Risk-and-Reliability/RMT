# Acceleration Factor Calculator
# Developed by Dr. Reuel Smith, 2021-2022

accelfactor <- function(params,lsm,Sacc,Suse){
  lsfcn<-lifestress.select(lsm)[[1]]
  if(length(Suse)==1 | lsm=="MultiStress") {
    AF<-lsfcn(params,Suse)/lsfcn(params,Sacc)
  }
  if(length(Suse)==2 & !(lsm=="MultiStress")) {
    AF<-lsfcn(params,Suse[1],Suse[2])/lsfcn(params,Sacc[1],Sacc[2])
  }
  return(AF)
}
