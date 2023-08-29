# Step-Stress Cumulative Damage Calculator
# Developed by Dr. Reuel Smith, 2021-2022

stepstress.damage <- function(data,stepstresstable){
  Ndat<-dim(data)[2]
  Ndat2<-dim(data)[1]
  Nss_tbl<-dim(stepstresstable)[2]
  Nss_tbl2<-dim(stepstresstable)[1]
  Nstress<-Ndat-2
  failcount<-rep(0,Nss_tbl2)

  # Initialize damage calculator
  dam<-rep(0,length(stepstresstable[,1]))
  if(Ndat==3){
    for(i2 in 1:Nss_tbl2){
      failcount[i2]<-sum(data[,2][which(data[,Ndat]==stepstresstable[i2,1])])
    }
  }else if(Ndat>=4){
    for(i2 in 1:Nss_tbl2){
      for(i3 in 1:Nstress){
        Smatch<-which(data[,i3+2]==stepstresstable[i2,i3])
        if(i3==1){
          Ssame<-Smatch
        } else{
          Ssame<-intersect(Smatch,Ssame)
        }
      }
      if(length(Ssame)==0){
        failcount[i2]<-0
      } else if(length(Ssame)>0){
        failcount[i2]<-sum(data[,2][Ssame])
      }
    }
  }
  dam<-cumsum(failcount)/Ndat2
  return(dam)
}
