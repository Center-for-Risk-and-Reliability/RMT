# Step-Stress Data Non-Cumulative Retabulator
# Developed by Dr. Reuel Smith, 2021-2022

stepstress.data <- function(data,stepstresstable) {
  library(dplyr)
  library(pracma)
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
  remnantstepstresstable<-integer(0)
  iremnant<-integer(0)

  # Check to see if rownames on stepstresstable exist.  They need to have names in order
  # for stepstress.LSQest to work properly
  if(is.null(rownames(stepstresstable))==TRUE){
    rownames(stepstresstable)<-num2str(c(1:Nss_tbl2),0)
  }

  # Correct or adjust the time if it is not already corrected
  TestTv<-rep(0,length(data[,Ndat]))
  CumTv<-rep(0,length(data[,Ndat]))
  cumtimes<-cumsum(stepstresstable[,Nss_tbl])
  if(Ndat==3){
    for (i2 in 1:length(TestTv)){
      TestTv[i2]<-stepstresstable[which(stepstresstable[,1]==data[i2,Ndat]),Nss_tbl]
      if(which(stepstresstable[,1]==data[i2,Ndat])-1>0){
        CumTv[i2]<-cumtimes[which(stepstresstable[,1]==data[i2,Ndat])-1]
      }
    }
  }else if(Ndat>=4){
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

  if(min(TestTv-data[,1])<0){
    adjT<-data[,1]-CumTv
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
  # Tally which times are in which step by stress
  Stepv <- rep(0,Ndat2)

  for (i2 in 1:length(stepstresstable[,1])) {
    if(Ndat==3){
      if(length(Sv)==0&&length(which(data[,Ndat]==stepstresstable[i2,1]))>0){
        Sv<-rep(stepstresstable[i2,1],length(which(data[,Ndat]>=stepstresstable[i2,1])))
        Cv<-c(data[,2][which(data[,Ndat]==stepstresstable[i2,1])],rep(0,length(which(data[,Ndat]>stepstresstable[i2,1]))))
        Tv<-c(adjT[which(data[,Ndat]==stepstresstable[i2,1])],rep(stepstresstable[i2,2],length(which(data[,Ndat]>stepstresstable[i2,1]))))
        iremnant<-i2
      } else if(length(Sv)>0&&length(which(data[,Ndat]==stepstresstable[i2,1]))>0){
        Sv<-c(Sv,rep(stepstresstable[i2,1],length(which(data[,Ndat]>=stepstresstable[i2,1]))))
        Cv<-c(Cv,data[,2][which(data[,Ndat]==stepstresstable[i2,1])],rep(0,length(which(data[,Ndat]>stepstresstable[i2,1]))))
        Tv<-c(Tv,adjT[which(data[,Ndat]==stepstresstable[i2,1])],rep(stepstresstable[i2,2],length(which(data[,Ndat]>stepstresstable[i2,1]))))
        iremnant<-c(iremnant,i2)
      }
      Stepv[which(data[,Ndat]== stepstresstable[i2,1])] <- i2
    } else if(Ndat>=4){
      for(i3 in 1:Nstress){
        Smatch<-which(data[,i3+2]==stepstresstable[i2,i3])
        if(i3==1){
          Ssame<-Smatch
        } else{
          Ssame<-intersect(Smatch,Ssame)
        }
      }
      if(length(Ssame)>0&&length(Sv[[1]])==0){
        for(i3 in 1:Nstress){
          Sv[[i3]]<-rep(stepstresstable[i2,i3],Ndat2)
        }
        Cv<-c(data[,2][Ssame],rep(0,Ndat2-last(Ssame)))
        Tv<-c(adjT[Ssame],rep(stepstresstable[i2,Nss_tbl],Ndat2-last(Ssame)))
        iremnant<-i2
      } else if(length(Ssame)>0&&length(Sv[[1]])>0){
        for(i3 in 1:Nstress){
          Sv[[i3]]<-c(Sv[[i3]],rep(stepstresstable[i2,i3],Ndat2-Ssame[1]+1))
        }
        Cv<-c(Cv,data[,2][Ssame],rep(0,Ndat2-last(Ssame)))
        Tv<-c(Tv,adjT[Ssame],rep(stepstresstable[i2,Nss_tbl],Ndat2-last(Ssame)))
        iremnant<-c(iremnant,i2)
      }
      Stepv[Ssame] <- i2
    }
  }
  if(Ndat==3){
    stress_v <- logical(0)
    for(i in 1:Nss_tbl2){
      if(i==1){
        stress_v <- rep(stepstresstable[i,1],Ndat2)
      }else{
        stress_v <- c(stress_v,rep(stepstresstable[i,1],Ndat2))
      }
    }
    stepstressdat <- matrix(c(Tv,Cv,Sv), nrow = length(Tv), ncol = Ndat, byrow = FALSE)
    updatedata0 <- matrix(c(rep(adjT,Nss_tbl2),rep(data[,2],Nss_tbl2),stress_v), nrow = Nss_tbl2*length(data[,2]), ncol = Ndat, byrow = FALSE)
  } else if (Ndat>=4){
    stress_v1 <- logical(0)
    stress_v2 <- logical(0)
    for(i in 1:Nss_tbl2){
      stress_v1 <- c(stress_v1,rep(stepstresstable[i,1],Ndat2))
      stress_v2 <- c(stress_v2,rep(stepstresstable[i,2],Ndat2))
    }
    stepstressdat <- matrix(c(Tv,Cv,unlist(Sv)), nrow = length(Tv), ncol = Ndat, byrow = FALSE)
    updatedata0 <- matrix(c(rep(adjT,Nss_tbl2),rep(data[,2],Nss_tbl2),stress_v1,stress_v2), nrow = Nss_tbl2*length(data[,2]), ncol = Ndat, byrow = FALSE)
  }
  return(list(stepstressdat,stepstresstable[iremnant,],updatedata0,Stepv))
}
