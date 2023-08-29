# Step-Stress Data Cumulative Retabulator
# Developed by Dr. Reuel Smith, 2021-2022

stepstress.data.cum <- function(data,stepstresstable) {
  library(dplyr)
  # Take the time data and the step stress table to set up your data
  # for step-stress analysis
  Ndat<-dim(data)[2]
  Ndat2<-dim(data)[1]
  Nss_tbl<-dim(stepstresstable)[2]
  Nss_tbl2<-dim(stepstresstable)[1]
  Nstress<-Ndat-2
  Sv<-integer(0)
  Cv<-integer(0)
  Tv<-integer(0)
  Stepv<-integer(0)
  Tendv<-integer(0)
  remnantstepstresstable<-integer(0)
  iremnant<-integer(0)

  # Correct or adjust the time if it is not already corrected
  TestTv<-rep(0,length(data[,Ndat]))
  CumTv<-rep(0,length(data[,Ndat]))
  cumtimes<-cumsum(stepstresstable[,Nss_tbl])
  if(Ndat==3){
    # CASE WHEN STRESS TYPE IS ONE
    for (i2 in 1:length(TestTv)){
      TestTv[i2]<-stepstresstable[which(stepstresstable[,1]==data[i2,Ndat]),Nss_tbl]
      if(which(stepstresstable[,1]==data[i2,Ndat])-1>0){
        CumTv[i2]<-cumtimes[which(stepstresstable[,1]==data[i2,Ndat])-1]
      }
    }
  }else if(Ndat>=4){
    # CASE WHEN STRESS TYPE IS TWO OR GREATER
    for (i2 in 1:Nss_tbl2){
      for(i3 in 1:Nstress){
        Smatch<-which(data[,i3+2]==stepstresstable[i2,i3])
        if(i3==1){
          Ssame<-Smatch
        } else{
          Ssame<-intersect(Smatch,Ssame)
        }
      }
      if(length(Ssame)>0){
        TestTv[Ssame]<-stepstresstable[i2,Nss_tbl]
        if(i2==1){
          CumTv[Ssame]<-0
        } else {
          CumTv[Ssame]<-cumtimes[i2-1]
        }
      }
    }
  }

  if(min(TestTv-data[,1])>=0){
    adjT<-data[,1]+CumTv
  } else{
    adjT<-data[,1]
  }
  # Identify if stress is single or multiple.  If multiple then set up a list
  # for each stress type.
  if(Ndat>=4){
    Slist<-vector("list",Nstress)
    Sv<-vector("list",Nstress)
    for(i2 in 1:Nstress){
      Slist[[i2]]<-data[,i2+2]
    }
  }

  for (i2 in 1:length(stepstresstable[,1])) {
    if(Ndat==3){
      # CASE WHEN STRESS TYPE IS ONE
      if(length(Sv)==0){
        Sv<-rep(stepstresstable[i2,1],length(which(data[,Ndat]>=stepstresstable[i2,1])))
        Cv<-c(data[,2][which(data[,Ndat]==stepstresstable[i2,1])],rep(0,length(which(data[,Ndat]>stepstresstable[i2,1]))))
        Tv<-c(adjT[which(data[,Ndat]==stepstresstable[i2,1])],rep(cumtimes[i2],length(which(data[,Ndat]>stepstresstable[i2,1]))))
        Stepv<-rep(i2,length(which(data[,Ndat]>=stepstresstable[i2,1])))
        if(i2==1){
          Tendv<-rep(0,length(which(data[,Ndat]>=stepstresstable[i2,1])))
        } else{
          Tendv<-rep(cumtimes[i2-1],length(which(data[,Ndat]>=stepstresstable[i2,1])))
        }
        iremnant<-i2
      } else if(length(Sv)>0){
        Sv<-c(Sv,rep(stepstresstable[i2,1],length(which(data[,Ndat]>=stepstresstable[i2,1]))))
        Cv<-c(Cv,data[,2][which(data[,Ndat]==stepstresstable[i2,1])],rep(0,length(which(data[,Ndat]>stepstresstable[i2,1]))))
        Tv<-c(Tv,adjT[which(data[,Ndat]==stepstresstable[i2,1])],rep(cumtimes[i2],length(which(data[,Ndat]>stepstresstable[i2,1]))))
        Stepv<-c(Stepv,rep(i2,length(which(data[,Ndat]>=stepstresstable[i2,1]))))
        if(i2==1){
          Tendv<-c(Tendv,rep(0,length(which(data[,Ndat]>=stepstresstable[i2,1]))))
        } else{
          Tendv<-c(Tendv,rep(cumtimes[i2-1],length(which(data[,Ndat]>=stepstresstable[i2,1]))))
        }
        iremnant<-c(iremnant,i2)
      }
    } else if(Ndat>=4){
      # CASE WHEN STRESS TYPE IS TWO OR GREATER
      Ssame<-integer(0)
      for(i3 in 1:Nstress){
        Smatch<-which(data[,i3+2]==stepstresstable[i2,i3])
        if(i3==1){
          Ssame<-Smatch
        } else{
          Ssame<-intersect(Smatch,Ssame)
        }
      }
      if(length(Ssame)>0){
        if(length(Sv[[1]])==0){
          for(i3 in 1:Nstress){
            Sv[[i3]]<-rep(stepstresstable[i2,i3],Ndat2)
          }
          Cv<-c(data[,2][Ssame],rep(0,Ndat2-last(Ssame)))
          Tv<-c(adjT[Ssame],rep(cumtimes[i2],Ndat2-last(Ssame)))
          Stepv<-rep(i2,Ndat2)
          if(i2==1){
            Tendv<-rep(0,Ndat2)
          } else{
            Tendv<-rep(cumtimes[i2-1],Ndat2)
          }
          iremnant<-i2
        } else if(length(Sv[[1]])>0){
          for(i3 in 1:Nstress){
            Sv[[i3]]<-c(Sv[[i3]],rep(stepstresstable[i2,i3],Ndat2-Ssame[1]+1))
          }
          Cv<-c(Cv,data[,2][Ssame],rep(0,Ndat2-last(Ssame)))
          Tv<-c(Tv,adjT[Ssame],rep(cumtimes[i2],Ndat2-last(Ssame)))
          Stepv<-c(Stepv,rep(i2,Ndat2-Ssame[1]+1))
          if(i2==1){
            Tendv<-c(Tendv,rep(0,Ndat2-Ssame[1]+1))
          } else{
            Tendv<-c(Tendv,rep(cumtimes[i2-1],Ndat2-Ssame[1]+1))
          }
          iremnant<-c(iremnant,i2)
        }
      } else if(length(Ssame)==0){
        if(length(Sv[[1]])==0){
          for(i3 in 1:Nstress){
            Sv[[i3]]<-rep(stepstresstable[i2,i3],Ndat2)
          }
          Cv<-rep(0,Ndat2)
          Tv<-rep(cumtimes[i2],Ndat2)
          Stepv<-rep(i2,Ndat2)
          if(i2==1){
            Tendv<-rep(0,Ndat2)
          } else{
            Tendv<-rep(cumtimes[i2-1],Ndat2)
          }
          iremnant<-i2
        } else if(length(Sv[[1]])>0){
          for(i3 in 1:Nstress){
            Sv[[i3]]<-c(Sv[[i3]],rep(stepstresstable[i2,i3],Ndat2))
          }
          Cv<-c(Cv,rep(0,Ndat2))
          Tv<-c(Tv,rep(cumtimes[i2],Ndat2))
          Stepv<-c(Stepv,rep(i2,Ndat2))
          if(i2==1){
            Tendv<-c(Tendv,rep(0,Ndat2))
          } else{
            Tendv<-c(Tendv,rep(cumtimes[i2-1],Ndat2))
          }
          iremnant<-c(iremnant,i2)
        }
      }
    }
  }
  if(Ndat==3){
    stepstressdat <- matrix(c(Tv,Cv,Sv), nrow = length(Tv), ncol = Ndat, byrow = FALSE)
  } else if (Ndat>=4){
    stepstressdat <- matrix(c(Tv,Cv,unlist(Sv)), nrow = length(Tv), ncol = Ndat, byrow = FALSE)
  }
  Tendiv<-Tendv[which(Cv==1)]
  Tendjv<-Tendv[which(Cv==0)]
  Stepiv<-Stepv[which(Cv==1)]
  Stepjv<-Stepv[which(Cv==0)]
  return(list(stepstressdat,Tendv,Tendiv,Tendjv,stepstresstable[iremnant,],Stepiv,Stepjv))
  #return(list(stepstressdat,Tendv,Tendiv,Tendjv,stepstresstable[iremnant,],Stepiv,Stepjv,Stepv))
}
