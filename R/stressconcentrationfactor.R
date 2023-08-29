# Stress Concentration Factor Calculator
# Developed by Reuel Smith, 2022

stress.concentration.factor <- function(dimensions,geometry,loadtype = 1){
  # Calculates stress concentration factor Kt based on Peterson's Stress Concentration Factors
  # Default for 'loadtype' will be tension (1).  Can also have bending (2) and traverse bending (3)
  # =========================================
  # Rectangular Bar or Plate options
  # =========================================
  # RB1. Rectangular Bar with a Semi-circle Edge Notch (input width W and radius r)
  # Chart 2.9 in Peterson's Stress Concentration Factors
  if(geometry == "rect_1semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    if(dimensions$r/dimensions$W >= 0.25){
      stop('r/W must be less than 0.25 for this geometry.')
    }
    x <- dimensions$r/dimensions$W
    Kt <- 3.065 - 8.871*x + 14.036*(x^2) - 7.219*(x^3)
  }
  # RB2. Rectangular Bar with Opposite Semi-circle Edge Notches (input width W and radius r)
  # Chart 2.3 in Peterson's Stress Concentration Factors
  if(geometry == "rect_2semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    if(dimensions$r/dimensions$W >= 0.25){
      x <- dimensions$r/(dimensions$W-2*dimensions$r)
      Kt <- 1.110*(x^-0.417)
    }
    if(dimensions$r/dimensions$W < 0.25){
      x <- (2*dimensions$r)/dimensions$W
      Kt <- 3.065 - 3.472*x + 1.009*(x^2) + 0.405*(x^3)
    }
  }
  # RB3. Rectangular Bar with a U-Notch (input width W, radius r, and U depth d)
  # Chart 2.10 in Peterson's Stress Concentration Factors
  if(geometry == "rect_1U_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    if(dimensions$d/dimensions$W >= 0.3 || dimensions$d/dimensions$r < 0.5 || dimensions$d/dimensions$r > 20) {
      stop('d/W must be less than 0.3 for this geometry and 0.5 < d/r < 20.')
    }
    if(dimensions$d/dimensions$r == 1){
      x <- dimensions$r/dimensions$W
      Kt <- 3.065 - 8.871*x + 14.036*(x^2) - 7.219*(x^3)
    }
    if(dimensions$d/dimensions$r != 1){
      if(dimensions$d/dimensions$r < 2){
        C1 <- 0.907 + 2.125*sqrt(dimensions$d/dimensions$r) + 0.023*(dimensions$d/dimensions$r)
        C2 <- 0.701 - 11.289*sqrt(dimensions$d/dimensions$r) + 1.708*(dimensions$d/dimensions$r)
        C3 <- -0.672 + 18.754*sqrt(dimensions$d/dimensions$r) - 4.046*(dimensions$d/dimensions$r)
        C4 <- 0.175 - 9.759*sqrt(dimensions$d/dimensions$r) + 2.365*(dimensions$d/dimensions$r)
      }
      if(dimensions$d/dimensions$r >= 2){
        C1 <- 0.953 + 2.136*sqrt(dimensions$d/dimensions$r) - 0.005*(dimensions$d/dimensions$r)
        C2 <- -3.255 - 6.281*sqrt(dimensions$d/dimensions$r) + 0.068*(dimensions$d/dimensions$r)
        C3 <- 8.203 + 6.893*sqrt(dimensions$d/dimensions$r) + 0.064*(dimensions$d/dimensions$r)
        C4 <- -4.851 - 2.793*sqrt(dimensions$d/dimensions$r) - 0.128*(dimensions$d/dimensions$r)
      }
      x <- dimensions$d/dimensions$W
      Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
    }
  }
  # RB4. Rectangular Bar with Opposite Edge U-Notches (input width W, radius r, and U depth d)
  # Chart 2.4 in Peterson's Stress Concentration Factors
  if(geometry == "rect_2U_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    if(dimensions$d/dimensions$W >= 0.3 || dimensions$d/dimensions$r < 0.5 || dimensions$d/dimensions$r > 20) {
      stop('d/W must be less than 0.3 for this geometry and 0.5 < d/r < 20.')
    }
    if(dimensions$d/dimensions$r == 1){
      x <- (2*dimensions$d)/dimensions$W
      Kt <- 3.065 - 3.472*x + 1.009*(x^2) + 0.405*(x^3)
    }
    if(dimensions$d/dimensions$r != 1){
      if(dimensions$d/dimensions$r < 2){
        C1 <- 0.955 + 2.169*sqrt(dimensions$d/dimensions$r) + 0.081*(dimensions$d/dimensions$r)
        C2 <- -1.557 - 4.046*sqrt(dimensions$d/dimensions$r) + 1.032*(dimensions$d/dimensions$r)
        C3 <- 4.013 + 0.424*sqrt(dimensions$d/dimensions$r) - 0.748*(dimensions$d/dimensions$r)
        C4 <- -2.461 + 1.538*sqrt(dimensions$d/dimensions$r) - 0.236*(dimensions$d/dimensions$r)
      }
      if(dimensions$d/dimensions$r >= 2){
        C1 <- 1.037 + 1.991*sqrt(dimensions$d/dimensions$r) + 0.002*(dimensions$d/dimensions$r)
        C2 <- -1.886 - 2.181*sqrt(dimensions$d/dimensions$r) - 0.048*(dimensions$d/dimensions$r)
        C3 <- 0.649 + 1.086*sqrt(dimensions$d/dimensions$r) + 0.142*(dimensions$d/dimensions$r)
        C4 <- 1.218 - 0.922*sqrt(dimensions$d/dimensions$r) - 0.086*(dimensions$d/dimensions$r)
      }
      x <- (2*dimensions$d)/dimensions$W
      Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
    }
  }
  # RB5. Rectangular Bar with a V-Notch (input width W, radius r, V depth d, and V angle alp)
  # Chart 2.7? in Peterson's Stress Concentration Factors

  # RB6. Rectangular Bar with Opposite Edge V-Notches (input width W, radius r,, V depth d, and V angle alp)
  # Chart 2.7 in Peterson's Stress Concentration Factors

  # RB7. Plate with a Circular Symmetric Hole (input width W, short diameter a, and long diameter d)
  # Chart 4.1 in Peterson's Stress Concentration Factors
  if(geometry == "rect_circ_symm_hole" && length(dimensions$W) == 1 && length(dimensions$d) == 1){
    if(dimensions$r/dimensions$W >= 0.25){
      stop('r/W must be less than 0.25 for this geometry.')
    }
    x <- 1 - dimensions$d/dimensions$W
    Kt <- 2 + 0.284*x - 0.6*(x^2) + 1.32*(x^3)
  }
  # RB8. Plate with a Elliptical Symmetric Hole (input width W and diameter d)
  # Chart 4.1 in Peterson's Stress Concentration Factors

  # RB9. Plate with a Slot Hole (input width W, slot distance d, and radius r)
  # Chart 4.1 in Peterson's Stress Concentration Factors

  # RB10. Plate with a Circular Offset Hole (input width W, edge to hole center distance a, and diameter d)
  # Chart 4.2 in Peterson's Stress Concentration Factors

  # RB11. Plate with a Deep Hyperbolic Notch (input edge to hyperbola distance d and thickness t)
  # Chart 4.2 in Peterson's Stress Concentration Factors

  # RB12. Plate with Two Deep Hyperbolic Notches (input edge to hyperbola distance d and thickness t)
  # Chart 2.1 in Peterson's Stress Concentration Factors

  # =========================================
  # Round or Shaft Options
  # =========================================
  # =========================================
  # Infinite/Semi-Infinite Plate Options
  # =========================================
  return(Kt)
}
