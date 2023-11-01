# Failure Time/Right-Censored Data sort (by Multiple-Stress)
# Developed by Dr. Reuel Smith, 2021-2022

sort.xircstressdata <- function(rawdat) {
  # Input should include stress and time data
  Nstress<-dim(rawdat)[2]-2

  if(length(rawdat[,2])==sum(rawdat[,2])) {
    xi <- rawdat[,1]
    if(Nstress==1){
      Sxi <- rawdat[,3:length(rawdat[1,])]
    } else {
      Sxi <- rawdat[,3:length(rawdat[1,])]
    }
    xitab<-list(xi,NULL,Sxi, NULL)
  }
  else {
    xi <- rawdat[,1][which(rawdat[,2] == 1)]
    rc <- rawdat[,1][which(rawdat[,2] == 0)]
    if(Nstress==1){
      Sxi <- rawdat[,3:length(rawdat[1,])][which(rawdat[,2] == 1)]
      Src <- rawdat[,3:length(rawdat[1,])][which(rawdat[,2] == 0)]
    } else {
      Sxilist <- vector("list",Nstress)
      Srclist <- vector("list",Nstress)
      for(i in 1:Nstress){
        Sxilist[[i]] <- rawdat[,i+2][which(rawdat[,2] == 1)]
        Srclist[[i]] <- rawdat[,i+2][which(rawdat[,2] == 0)]
      }
      Sxi <- matrix(unlist(Sxilist),nrow = length(Sxilist[[1]]), ncol=Nstress,byrow = FALSE)
      Src <- matrix(unlist(Srclist),nrow = length(Srclist[[1]]), ncol=Nstress,byrow = FALSE)
    }
    xitab<-list(xi,rc,Sxi,Src)
  }
  return(xitab)
}
