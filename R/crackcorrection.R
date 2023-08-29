# Crack Correction Factor Calculator
# Developed by Reuel Smith, 2022

crack.correction <- function(dimensions,crackgeometry){
  # Computes correction factor or f(g) for a specimen of a given geometry
  if(crackgeometry == "center_2a" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    fg <-sqrt(1/cos((pi*dimensions$r)/dimensions$W))
  }

  if(crackgeometry == "edge_single_1a" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    x <- dimensions$r/dimensions$W
    fg <- 1.12 - 0.231*x + 10.55*(x^2) - 21.72*(x^3) + 30.39*(x^4)
  }

  if(crackgeometry == "edge_double_2a" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    x <- dimensions$r/(0.5*dimensions$W)
    fg <- 1.12 + 0.203*x - 1.197*(x^2) + 1.93*(x^3)
  }
  if(crackgeometry == "edge_semi_circle_thick_body_1a" || crackgeometry == "corner_circle_thick_body_1a"){
    fg <- 1.12
  }

  return(fg)
}
