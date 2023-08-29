# Data Sort by Stress
# Developed by Dr. Reuel Smith, 2021-2022

checkstress <- function(dat) {
  # Checks data block for third column or higher
  # Parses the data into blocks for easy probability plotting
  # Outputs the data as a list
  colsize<-dim(dat)[2]
  if(colsize == 2){
    # No stress case
    dataout<-dat
  }
  if(colsize == 3){
    # Case where one stress type exists
    stressset<-dat[,3][!duplicated(dat[,3])]
    dataout<-vector(mode = "list", length = length(stressset))
    for(i in 1:length(stressset)) {
      dataout[[i]]<-dat[which(dat[,3] == stressset[i]), ]
    }
  }
  if(colsize >= 4){
    # Case where multiple stresses exist
    stressset<-dat[,3:colsize][!duplicated(dat[,3:colsize]), ]
    dataout<-vector(mode = "list", length = length(stressset[,1]))
    datcheck<-rowMeans(dat[,3:colsize])
    stresscheck<-rowMeans(stressset)
    for(i in 1:length(stressset[,1])) {
      dataout[[i]]<-dat[which(datcheck == stresscheck[i]), ]
    }
  }
  return(dataout)
}
