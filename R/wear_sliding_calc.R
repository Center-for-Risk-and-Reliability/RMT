# Sliding Wear Calculator
# Developed by Reuel Smith, 2023

wear.sliding <- function(data,matproperties = list(k = 0, mat1 = "MildSteel", mat2 = "MildSteel"),units = 1){
  library(pracma)

  kdatabase = data.frame(material1 = c("Zinc", "MildSteel", "Platinum", "Copper", "StainlessSteel", "Silver", "Gold",
                                       "Copper", "Platinum", "Platinum", "Brass70_30", "Bakelite", "MildSteel"),
                         material2 = c("Zinc", "MildSteel", "Platinum", "Copper", "StainlessSteel", "Silver", "Gold",
                                       "MildSteel", "MildSteel", "Silver", "Steel4140", "Bakelite", "Copper"),
                         k = c(53.3333333e-3, 15e-3, 13e-3, 10.666666667e-3, 7e-3, 4e-3, 0.1,
                               0.5e-3, 0.5e-3, 0.333333333333e-3, 4.3e-4, 0.02e-3, 1.6666667e-4))

  # Check units and set up axis labels
  if(units == 1){
    # For metric units, input must be: wear length (mm), stress (MPa), force (N),
    velocitylabel <- " mm\U00B3/sec"
    wearlabel <- " mm\U00B3"
    velcorr <- 1000^2
  }
  if(units == 2){
    # For English units, input must be: wear length (inches), stress (ksi), force (kips)
    velocitylabel <- " in\U00B3/sec"
    wearlabel <- " in\U00B3"
    velcorr <-  1
  }

  # Initialize wear coefficient k from mat1 (material 1) on mat2 (material 2)
  if (matproperties$k == 0 || length(matproperties$k) == 0){
    k <- kdatabase$k[intersect(which(kdatabase$material1 == matproperties$mat1),which(kdatabase$material2 == matproperties$mat2))]
  } else{
    k <- matproperties$k
  }

  # Data Input Type #1: k, load (N), sliding distance (L), and minimum hardness (H)
  # RETURN WEAR
  if(length(data$P) > 0 && length(data$V) == 0 && length(matproperties$H) > 0 && length(data$A) == 0 && length(data$h) == 0 && length(data$L) > 0){
    W <- k*((data$P*data$L)/min(matproperties$H))
    output <- list(wear = W)
    # Produce some output text that summarizes the results
    cat(c("Sliding/Adhesive Wear is ",W,wearlabel,"\n\n"))
  }

  # Data Input Type #2: Area (A) and wear depth (h)
  # RETURN WEAR
  if(length(data$P) == 0 && length(data$V) == 0 && length(matproperties$H) == 0 && length(data$A) > 0 && length(data$h) > 0  && length(data$L) == 0){
    W <- data$A*data$h
    output <- list(wear = W)
    # Produce some output text that summarizes the results
    cat(c("Sliding/Adhesive Wear is ",W,wearlabel,"\n\n"))
  }

  # Data Input Type #3: k, load (N), Cutting speed (V), and minimum hardness (H)
  # RETURN WEAR VELOCITY
  if(length(data$P) > 0 && length(data$V) > 0 && length(matproperties$H) > 0 && length(data$A) == 0 && length(data$h) == 0  && length(data$L) == 0){
    W_vel <- k*((data$P*data$V)/min(matproperties$H))*velcorr
    output <- list(wear_Velocity = W_vel)
    # Produce some output text that summarizes the results
    cat(c("Sliding/Adhesive Wear Velocity is ",W_vel,velocitylabel,"\n\n"))
  }

  # Data Input Type #4: k, load (N), Cutting speed (V), minimum hardness (H), and wear depth (h)
  # RETURN WEAR, WEAR VELOCITY, and TIME TO REMOVE wear h
  if(length(data$P) > 0 && length(data$V) > 0 && length(matproperties$H) > 0 && length(data$A) > 0 && length(data$h) > 0  && length(data$L) == 0){
    W_vel <- k*((data$P*data$V)/min(matproperties$H))*velcorr
    W <- data$A*data$h
    time <- W/W_vel
    output <- list(wear = W, wear_Velocity = W_vel, time_to_wear = time)
    # Produce some output text that summarizes the results
    cat(c("Sliding/Adhesive Wear is ",W,wearlabel,"\nSliding/Adhesive Wear Velocity is ",W_vel,velocitylabel,"\nElapsed time to remove Sliding/Adhesive Wear is ",time," sec\n\n"))
  }

  return(output)
}
