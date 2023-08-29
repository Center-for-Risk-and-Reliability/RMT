# Accelerated Degradation Testing Rank System
# Developed by Dr. Reuel Smith, 2021-2022

adt.rank <- function(data){
  # Short function that ranks ADT data based on appropriateness of destruction model

  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  R2_1<-adt.full.LSQ(data,"Linear",1)[[1]][,dim(adt.full.LSQ(data,"Linear",1)[[1]])[2]]
  R2_2<-adt.full.LSQ(data,"Exponential",1)[[1]][,dim(adt.full.LSQ(data,"Exponential",1)[[1]])[2]]
  R2_3<-adt.full.LSQ(data,"SquareRoot",1)[[1]][,dim(adt.full.LSQ(data,"SquareRoot",1)[[1]])[2]]
  R2_4<-adt.full.LSQ(data,"Power",1)[[1]][,dim(adt.full.LSQ(data,"Power",1)[[1]])[2]]
  R2_5<-adt.full.LSQ(data,"Logarithmic",1)[[1]][,dim(adt.full.LSQ(data,"Logarithmic",1)[[1]])[2]]
  # R2_6<-adt.full.LSQ(data,"Gompertz",1)[[1]][,dim(adt.full.LSQ(data,"Gompertz",1)[[1]])[2]]
  R2_7<-adt.full.LSQ(data,"LloydLipow",1)[[1]][,dim(adt.full.LSQ(data,"LloydLipow",1)[[1]])[2]]

  if(min(data[,2])>0 && max(data[,2])<1){
    R2_8<-adt.full.LSQ(data,"Mitsuom",0.5)[[1]][,dim(adt.full.LSQ(data,"Mitsuom",0.5)[[1]])[2]]

    # modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Gompertz","Lloyd-Lipow","Mitsuom")
    # R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7,R2_8), nrow = 8, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    # R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7,R2_8), nrow = 8, ncol = length(unitnames), byrow = TRUE)
    modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Lloyd-Lipow","Mitsuom")
    R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7,R2_8), nrow = 7, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7,R2_8), nrow = 7, ncol = length(unitnames), byrow = TRUE)
  } else{
    # modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Gompertz","Lloyd-Lipow")
    # R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7), nrow = 7, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    # R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7), nrow = 7, ncol = length(unitnames), byrow = TRUE)
    modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Lloyd-Lipow")
    R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7), nrow = 6, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7), nrow = 6, ncol = length(unitnames), byrow = TRUE)
  }

  # Rank and average rank
  for(i in 1:length(unitnames)){
    if(i==1){
      rankset<-dim(R2set)[1] + 1 - rank(R2setblank[,i], ties.method = "first")
    } else{
      rankset<-c(rankset,dim(R2set)[1] + 1 - rank(R2setblank[,i], ties.method = "first"))
    }
  }
  for(i in 1:length(modelnames)){
    avgR2 <- mean(R2set[i,])
    if(i==1){
      avgR2set <- avgR2
    } else{
      avgR2set <- c(avgR2set,avgR2)
    }
  }
  avgrankset<- dim(R2set)[1] + 1 - rank(avgR2set, ties.method = "first")
  avgrankset<- matrix(avgrankset,nrow =dim(R2set)[1], ncol = 1, byrow = FALSE, dimnames = list(modelnames,"Average Rank"))
  rankset <- matrix(rankset,nrow = dim(R2set)[1], ncol = dim(R2set)[2], byrow = FALSE, dimnames = list(modelnames,unitnames))

  return(list(R2set,rankset,avgrankset))
}
