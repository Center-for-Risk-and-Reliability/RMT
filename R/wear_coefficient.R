# Sliding Wear Coefficient Tool
# Developed by Reuel Smith, 2022

wear.coefficient <- function(mat1,mat2){
  # Selects wear coefficient from mat1 (material 1) on mat2 (material 2) contact list.
  # List includes: Zinc, mild steel, platinum, copper, stainless steel, silver

  # Zinc on Zinc
  if(mat1 == "Zinc" && mat2 == "Zinc"){
    k <- 53.3333333e-3
  }

  # Mild steel on Mild steel
  if(mat1 == "MildSteel" && mat2 == "MildSteel"){
    k <- 15e-3
  }

  # Platinum on Platinum
  if(mat1 == "Platinum" && mat2 == "Platinum"){
    k <- 13e-3
  }

  # Copper on Copper
  if(mat1 == "Copper" && mat2 == "Copper"){
    k <- 10.666666667e-3
  }
  # Stainless Steel on Stainless Steel
  if(mat1 == "StainlessSteel" && mat2 == "StainlessSteel"){
    k <- 7e-3
  }

  # Silver on Silver
  if(mat1 == "Silver" && mat2 == "Silver"){
    k <- 4e-3
  }

  # Gold on Gold (states that the range is from 0.1 to 1.  Will need to verify this.)
  if(mat1 == "Gold" && mat2 == "Gold"){
    k <- 0.1
  }

  # Copper on Mild Steel
  if(mat1 == "Copper" && mat2 == "MildSteel"){
    k <- 0.5e-3
  }

  # Platinum on Mild Steel
  if(mat1 == "Platinum" && mat2 == "MildSteel"){
    k <- 0.5e-3
  }

  # Platinum on Silver
  if(mat1 == "Platinum" && mat2 == "Silver"){
    k <- 0.333333333333e-3
  }

  # Steel 4140 on 70-30 Brass
  if(mat1 == "Brass70_30" && mat2 == "Steel4140"){
    k <- 4.3e-4
  }

  return(k)
}
