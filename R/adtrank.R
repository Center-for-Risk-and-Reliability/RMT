# Accelerated Degradation Testing Rank System
# Developed by Dr. Reuel Smith, 2021-2025

adt.rank <- function(data){
  # Short function that ranks ADT data based on appropriateness of destruction model
  # ===============================================================
  # UPDATE: RCS 02212024
  # New adt.rank to use degradationlife.LSQest
  # UPDATE: RCS 02222024
  # Now determines fitness by measure of Sum of squares error (SSE)
  # UPDATE: RCS 01132026
  # Updated to suit new output system of degradationlife.LSQest.  Also sets endurance limit to the maximum
  # (or minimum) degradation
  # ===============================================================
  # Pulls the unit designations from column 3 of the input data
  unitnames <- unique(data[,3])

  # Find slope of the first unit to get general direction of degradation (increasing or decreasing)
  xfit <- data[which(data[,3]==unitnames[1]),1]
  yfit <- data[which(data[,3]==unitnames[1]),2]
  params  <- lm(yfit ~ poly(xfit, 1, raw=TRUE))
  slope <- summary(params)$coefficients[2,1]

  if(slope > 0){ # Set endurance limit to max degradation if increasing
    D0 <- max(data[,2])
  }
  if(slope < 0){ # Set endurance limit to min degradation if decreasing
    D0 <- min(data[,2])
  }


  R2_1<-degradationlife.LSQest(data,"Linear",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Linear",D0=D0)[[1]])[2]]
  R2_2<-degradationlife.LSQest(data,"Exponential",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Exponential",D0=D0)[[1]])[2]]
  R2_3<-degradationlife.LSQest(data,"SquareRoot",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"SquareRoot",D0=D0)[[1]])[2]]
  R2_4<-degradationlife.LSQest(data,"Power",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Power",D0=D0)[[1]])[2]]
  R2_5<-degradationlife.LSQest(data,"Logarithmic",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Logarithmic",D0=D0)[[1]])[2]]
  # R2_6<-degradationlife.LSQest(data,"Gompertz",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Gompertz",D0=D0)[[1]])[2]]
  R2_7<-degradationlife.LSQest(data,"LloydLipow",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"LloydLipow",D0=D0)[[1]])[2]]
  # R2_8<-degradationlife.LSQest(data,"Mitsuom",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Mitsuom",D0=D0)[[1]])[2]]
  R2_9<-degradationlife.LSQest(data,"SquareRoot2",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"SquareRoot2",D0=D0)[[1]])[2]]


  if(min(data[,2])>0 && max(data[,2])<1){
    R2_8<-degradationlife.LSQest(data,"Mitsuom",D0=D0)[[1]][,dim(degradationlife.LSQest(data,"Mitsuom",D0=D0)[[1]])[2]]

    # modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Gompertz","Lloyd-Lipow","Mitsuom")
    # R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7,R2_8), nrow = 8, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    # R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7,R2_8), nrow = 8, ncol = length(unitnames), byrow = TRUE)
    modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Lloyd-Lipow","Mitsuom","Square-Root2")
    R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7,R2_8,R2_9), nrow = 8, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7,R2_8,R2_9), nrow = 8, ncol = length(unitnames), byrow = TRUE)
  } else{
    # modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Gompertz","Lloyd-Lipow")
    # R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7), nrow = 7, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    # R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6,R2_7), nrow = 7, ncol = length(unitnames), byrow = TRUE)
    modelnames <- c("Linear","Exponential","Square-Root","Power","Logarithmic","Lloyd-Lipow","Square-Root2")
    R2set <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7,R2_9), nrow = 7, ncol = length(unitnames), byrow = TRUE, dimnames = list(modelnames,unitnames))
    R2setblank <- matrix(c(R2_1,R2_2,R2_3,R2_4,R2_5,R2_7,R2_9), nrow = 7, ncol = length(unitnames), byrow = TRUE)
  }

  # Rank and average rank
  for(i in 1:length(unitnames)){
    if(i==1){
      rankset<-rank(R2setblank[,i], ties.method = "first")
      # rankset<-dim(R2set)[1] + 1 - rank(R2setblank[,i], ties.method = "first")
    } else{
      rankset<-c(rankset, rank(R2setblank[,i], ties.method = "first"))
      # rankset<-c(rankset,dim(R2set)[1] + 1 - rank(R2setblank[,i], ties.method = "first"))
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
  avgrankset<- rank(avgR2set, ties.method = "first")
  # avgrankset<- dim(R2set)[1] + 1 - rank(avgR2set, ties.method = "first")
  avgrankset<- matrix(avgrankset,nrow =dim(R2set)[1], ncol = 1, byrow = FALSE, dimnames = list(modelnames,"Average Rank"))
  rankset <- matrix(rankset,nrow = dim(R2set)[1], ncol = dim(R2set)[2], byrow = FALSE, dimnames = list(modelnames,unitnames))

  return(list(SSE.by.unit = R2set,adt.rank.by.unit = rankset,adt.rank.average = avgrankset))
}
