# Data Count by Stress
# Developed by Dr. Reuel Smith, 2020-2022

checkdatacount <- function(databystress){
  # Load 'pracma' library to use 'size' function
  library(pracma)
  # Check if the data is multi-stress or not to continue.
  if(is.null(size(databystress))==FALSE){
    # Setup vector of stress cases where data is greater than or equal to 2
    datcountGTE2<-rep(0,size(databystress)[2])

    # Check for Censored data too.  If no censor cases are 1, then set as zero
    censcount<-rep(0,size(databystress)[2])
    for(i in 1:size(databystress)[2]){
      # Check row count of each list entry.  Need to alternate between use of 'DIM' and
      # 'SIZE' because row name may affect what is pulled depending on the data
      if(is.null(rownames(databystress[[i]])) == TRUE){
        # Case where there are no row names in list data entry
        if(size(databystress[[i]])[1] >= 2){
          datcountGTE2[i]<-1
          if(sum(databystress[[i]][,2]) >= 1){
            censcount[i]<-1
          }
        } else {
          if(sum(databystress[[i]][2])>=1){
            censcount[i]<-1
          }
        }
      } else{
        # Case where there ARE row names in list data entry
        if(dim(databystress[[i]])[1] >= 2){
          datcountGTE2[i]<-1
          if(sum(databystress[[i]][,2]) >= 1){
            censcount[i]<-1
          }
        } else {
          if(sum(databystress[[i]][2])>=1){
            censcount[i]<-1
          }
        }
      }
    }
  } else{
    datcountGTE2<-c(1)
    censcount<-c(1)
  }
  return(list(datcountGTE2,censcount))
}
