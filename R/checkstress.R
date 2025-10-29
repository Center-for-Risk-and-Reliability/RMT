# Data Sort by Stress
# Developed by Dr. Reuel Smith, 2021-2025

checkstress <- function(dat) {
  # Checks data block for third column or higher (fourth block and higher for ADT data)
  # Parses the data into blocks for easy probability plotting
  # Outputs the data as a list

  # RCS04022025 - Modify to be able to read ADT data as well as ALT data.

  colsize<-dim(dat)[2] # Pulls the number of columns of data block

  if(colsize == 2){
    # No stress case
    dataout<-dat
  }

  # Check to see if data block is for ALT or ADT.  ADT will have a non-integer descriptor in its unit column (column 3)
  if (colsize >= 3 && is.character(data.frame(dat)[1,3]) == FALSE){ # It is ALT if column 3 is a double
    COL_0 <- 3
    COL_1 <- 4
  }
  if (colsize >= 4 && is.character(data.frame(dat)[1,3]) == TRUE){ # It is ADT if column 3 is not a double
    COL_0 <- 4
    COL_1 <- 5
  }

  if(colsize == COL_0){
    # Case where one stress type exists
    stressset<-dat[,COL_0][!duplicated(dat[,COL_0])]
    dataout<-vector(mode = "list", length = length(stressset))
    for(i in 1:length(stressset)) {
      dataout[[i]]<-dat[which(dat[,COL_0] == stressset[i]), ]
    }
  }
  if(colsize >= COL_1){
    # Case where multiple stresses exist
    stressset<-dat[,COL_0:colsize][!duplicated(dat[,COL_0:colsize]), ]
    dataout<-vector(mode = "list", length = length(stressset[,1]))
    datcheck<-rowMeans(dat[,COL_0:colsize])
    stresscheck<-rowMeans(stressset)
    for(i in 1:length(stressset[,1])) {
      dataout[[i]]<-dat[which(datcheck == stresscheck[i]), ]
    }
  }
  return(dataout)
}
