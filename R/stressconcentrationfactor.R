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
  # Chart 2.9 in Peterson's Stress Concentration Factors
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
  # # =========================================
  # # Camille's Contributions to SCF (test now Sort later)
  # # =========================================
  # # RB6. Rectangular Bar with Opposite Edge V-Notches (input width W, radius r, V depth d, and V angle alp)
  # # Chart 2.7 in Peterson's Stress Concentration Factors
  # if(geometry == "rect_2V_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1 && length(dimensions$alp) == 1){
  #   if(dimensions$alp > 150  || !(setequal((2*dimensions$d/dimensions$w),0.4) || !(setequal((2*dimensions$d/dimensions$w),2/3))) || dimensions$d/dimensions$r <= 1 || dimensions$d/dimensions$r > 50){
  #     stop('alpha < 150, 2d/w = 0.4 or 0.667, and 1 < d/r < 50 for this geometry.')
  #   }
  #   if(dimensions$d/dimensions$r < 2){
  #     C1 <- 0.955 + 2.169*sqrt(dimensions$d/dimensions$r) + 0.081*(dimensions$d/dimensions$r)
  #     C2 <- -1.557 - 4.046*sqrt(dimensions$d/dimensions$r) + 1.032*(dimensions$d/dimensions$r)
  #     C3 <- 4.013 + 0.424*sqrt(dimensions$d/dimensions$r) - 0.748*(dimensions$d/dimensions$r)
  #     C4 <- -2.461 + 1.538*sqrt(dimensions$d/dimensions$r) - 0.236*(dimensions$d/dimensions$r)
  #   }
  #   if(dimensions$d/dimensions$r >= 2){
  #     C1 <- 1.037 + 1.991*sqrt(dimensions$d/dimensions$r) + 0.002*(dimensions$d/dimensions$r)
  #     C2 <- -1.886 - 2.181*sqrt(dimensions$d/dimensions$r) - 0.048*(dimensions$d/dimensions$r)
  #     C3 <- 0.649 + 1.086*sqrt(dimensions$d/dimensions$r) + 0.142*(dimensions$d/dimensions$r)
  #     C4 <- 1.218 - 0.922*sqrt(dimensions$d/dimensions$r) - 0.086*(dimensions$d/dimensions$r)
  #   }
  #   x <- (2*dimensions$d)/dimensions$W
  #   Ktu <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  #
  #   if(2*dimensions$d/dimensions$w == 0.4){
  #     if(Ktu > 3.5 || Ktu < 1.6){
  #       stop('1.6 < Ktu < 3.5 for these parameters.')
  #     }
  #     if(dimensions$alp < 90){
  #       Kt <- Ktu
  #     }
  #
  #     CC1 <- 5.294 - 0.1225*dimensions$alp + 0.000523*(dimensions$alp^2)
  #     CC2 <--5.0002 + 0.1171*dimensions$alp - 0.000434*(dimensions$alp^2)
  #     CC3 <-1.423 - 0.00197*dimensions$alp - 0.000004*(dimensions$alp^2)
  #
  #     Kt <- CC1 + CC2*(Ktu^0.5) + CC3*Ktu
  #   }
  #
  #   if(2*dimensions$d/dimensions$w == 2/3){
  #     if(Ktu > 2.8 || Ktu < 1.6){
  #       stop('1.6 < Ktu < 2.8 for these parameters.')
  #     }
  #     if(dimensions$alp < 60){
  #       Kt <- Ktu
  #     }
  #     CC1 <- -10.01 + 0.1534*dimensions$alp - 0.000647*(dimensions$alp^2)
  #     CC2 <- 13.60 - 0.2140*dimensions$alp + 0.000973*(dimensions$alp^2)
  #     CC3 <- -3.781 + 0.07875*dimensions$alp - 0.000392*(dimensions$alp^2)
  #
  #     Kt <- CC1 + CC2*(Ktu^0.5) + CC3*Ktu
  #   }
  # }
  #
  # # RB7. Rectangular Bar with Infinite Rows of Semicircular Notches (input width w, spacing b, radius r)
  # # Chart 2.12 in Peterson's Stress Concentration Factors
  # if(geometry == "rect_inf_semicirc_edge" && length(dimensions$w) == 1 && length(dimensions$b) == 1 && length(dimensions$r) == 1){
  #   if(dimensions$r/dimensions$w > 0.4 || dimensions$r/dimensions$b > 1){
  #     stop('r/w must be less than 0.4 for this geometry and r/b must be < 1.')
  #   }
  #   x <- (dimensions$r/dimensions$w)
  #   C1 <- 3.1055 - 3.2487*x + 0.8522*(x^2)
  #   C2 <- - 1.437 + 10.5053*x - 8.7547*(x^2) + 19.6273*(x^3)
  #   C3 <- -1.6753 - 14.0851*x + 43.675*(x^2)
  #   C4 <- 1.7207 + 5.7974*x - 27.7453*(x^2) + 6.0444*(x^3)
  #
  #   y <- (dimensions$r/dimensions$b)
  #   Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  # }
  #
  # # =========================================
  # # Infinite/Semi-Infinite Plate Options
  # # =========================================
  #
  # # P1. Plate with a Circular Symmetric Hole (input width W, diameter d)
  # # Chart 4.1 in Peterson's Stress Concentration Factors
  # if(geometry == "rect_circ_symm_hole" && length(dimensions$W) == 1 && length(dimensions$d) == 1){
  #   if(dimensions$d*0.5/dimensions$W >= 0.25){
  #     stop('r/W must be less than 0.25 for this geometry.')
  #   }
  #   x <- 1 - dimensions$d/dimensions$W
  #   Kt <- 2 + 0.284*x - 0.6*(x^2) + 1.32*(x^3)
  # }
  #
  # # P2. Plate with an Elliptical Symmetric Hole (input width W, major axis a, and minor axis b)
  # # Chart 4.51 in Peterson's Stress Concentration Factors
  # if(geometry == "plate_ellip_symm_hole" && length(dimensions$w) == 1 && length(dimensions$a) == 1 && length(dimensions$b) == 1){
  #   if(dimensions$a/dimensions$w > 0.5 || dimensions$a/dimensions$b < 1 || dimensions$a/dimensions$b > 8){
  #     stop('a/w must be less than 0.5 for this geometry and 1 < a/b < 8.')
  #   }
  #   x <- dimensions$a/dimensions$b
  #   C1 <- 1.109 - 0.188*(x^0.5) + 2.086*x
  #   C2 <- -0.486 + 0.213*(x^0.5) - 2.588*x
  #   C3 <- 3.816 - 5.51*(x^0.5) + 4.638*x
  #   C4 <- -2.438 + 5.485*(x^0.5) -4.126*x
  #
  #   y <- 2*dimensions$a/dimensions$w
  #
  #   Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  # }
  #
  # # P3. Plate with a Slot Hole (input width W, slot distance d, and radius r)
  # # Chart 4.1 in Peterson's Stress Concentration Factors
  #
  # # P4. Plate with a Circular Offset Hole (input width W, edge to hole center distance a, and radius r)
  # # Chart 4.3 in Peterson's Stress Concentration Factors
  # if(geometry == "plate_circ_offset_hole" && length(dimensions$w) == 1 && length(dimensions$a) == 1 && length(dimensions$r) == 1){
  #   if(dimensions$r/dimensions$a > 0.5 || (dimensions$w-dimensions$a)/dimensions$a < 1){
  #     stop('r/a must be less than 0.5 for this geometry and (w-a)/a < 1.')
  #   }
  #   x <- (dimensions$a/(dimensions$w-dimensions$a))
  #   C1 <- 2.989-0.0064*x
  #   C2 <- -2.872 + 0.095*x
  #   C3 <- 2.348 + 0.196*x
  #
  #   y <- (dimensions$r/dimensions$a)
  #   Kt <- C1 + C2*y + C3*(y^2)
  # }
  #
  # # P5. Plate with a Deep Hyperbolic Notch (input edge to hyperbola distance d and thickness t)
  # # Chart 4.2 in Peterson's Stress Concentration Factors
  #
  # # P6. Plate with Two Deep Hyperbolic Notches (input edge to hyperbola distance d and thickness t)
  # # Chart 2.1 in Peterson's Stress Concentration Factors
  #
  # # =========================================
  # # Shaft Options
  # # =========================================
  #
  # #S1. Shaft with a Shallow U-Shaped Groove, Tension (input diameter D, radius of groove r, depth of groove t)
  # #Chart 2.21 in Peterson's Stress Concentration Factors
  # if(geometry == "shaft_wideU_groove_tens" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$t) == 1){
  #   d <- dimensions$D-(2*dimensions$t)
  #   if(dimensions$D/d > 1.1 || dimensions$D/d < 1.005 || dimensions$r/d > 1 || dimensions$r/d < 0.3){
  #     stop('1.005 < D/(D-2t) < 1.1 for this geometry and 0.3 < r/(D-2t) < 1.')
  #   }
  #   x <- dimensions$D/d
  #   C1 <- -81.39 + 153.1*x - 70.49*(x^2)
  #   C2 <- 119.64 - 221.81*x + 101.93*(x^2)
  #   C3 <- -57.88 + 107.33*x - 49.34*(x^2)
  #
  #   y <- dimensions$r/d
  #   Kt <- C1 + C2*y + C3*(y^2)
  # }
  #
  # #S2. Shaft with a Shallow U-Shaped Groove, Torsion (input diameter D, radius of groove r, depth of groove t)
  # #Chart 2.50 in Peterson's Stress Concentration Factors
  # if(geometry == "shaft_wideU_groove_tors" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$t) == 1){
  #   d <- dimensions$D-(2*dimensions$t)
  #   if(dimensions$D/d > 1.1 || dimensions$D/d < 1.005 || dimensions$r/d > 1 || dimensions$r/d < 0.3){
  #     stop('1.005 < D/(D-2t) < 1.1 for this geometry and 0.3 < r/(D-2t) < 1.')
  #   }
  #   x <- dimensions$D/d
  #   C1 <- -35.16 + 67.57*x - 31.28*(x^2)
  #   C2 <- 79.13 - 148.37*x + 69.09*(x^2)
  #   C3 <- -50.34 + 94.67*x - 44.26*(x^2)
  #
  #   y <- dimensions$r/d
  #   Kt <- C1 + C2*y + C3*(y^2)
  # }
  #
  # #S3. Shaft with a U-Shaped or Semicircular Groove, Torsion (input diameter D, radius of groove r, depth of groove t)
  # #Chart 2.47 in Peterson's Stress Concentration Factors
  # if(geometry == "shaft_regU_groove_tors" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$t) == 1){
  #   x <- dimensions$t/dimensions$r
  #   if(x < 0.25 || x > 50){
  #     stop('0.25 < t/r < 50 for this geometry.')
  #   }
  #   y <- 2*dimensions$t/dimensions$D
  #
  #   if(x == 1){
  #     C1 <- 2
  #     C2 <- -3.555
  #     C3 <- 4.898
  #     C4 <- -2.365
  #   }
  #   if(x !=1 || x < 2){
  #     C1 <- 0.966 + 1.056*(x^0.5) - 0.022*x
  #     C2 <- -0.192 - 4.037*(x^0.5) + 0.674*x
  #     C3 <- 0.808 + 5.321*(x^0.5) - 1.231*x
  #     C4 <- -0.567 - 2.364*(x^0.5) + 0.566*x
  #   }
  #   if(x >= 2){
  #     C1 <- 1.089 + 0.924*(x^0.5) + 0.018*x
  #     C2 <- -1.504 - 2.141*(x^0.5) - 0.047*x
  #     C3 <- 2.486 + 2.289*(x^0.5) + 0.091*x
  #     C4 <- 1.056 - 1.104*(x^0.5) - 0.059*x
  #   }
  #
  #   Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  # }
  #
  # #S4. Shaft with a V-Shaped Groove, Torsion (input diameter D, radius r, V depth d, and V angle alp)
  # #Chart 2.51 in Peterson's Stress Concentration Factors
  # if(geometry == "shaft_V_groove_tors" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1 && length(dimensions$alp) == 1){
  #   dia <- dimensions$D - (2*dimensions$d)
  #
  #   if(dimensions$alp > 125 || (dimensions$alp > 90 && dimensions$r/dia > 0.01)){
  #     stop('alpha < 125 for this geometry, and for alpha > 90, r/(D-2d) < 0.01.')
  #   }
  #
  #   x <- dimensions$d/dimensions$r
  #   if(x < 0.25 || x > 50){
  #     stop('0.25 < d/r < 50 for this geometry.')
  #   }
  #   y <- 2*dimensions$d/dimensions$D
  #
  #   if(x == 1){
  #     C1 <- 2
  #     C2 <- -3.555
  #     C3 <- 4.898
  #     C4 <- -2.365
  #   }
  #   if(x !=1 || x < 2){
  #     C1 <- 0.966 + 1.056*(x^0.5) - 0.022*x
  #     C2 <- -0.192 - 4.037*(x^0.5) + 0.674*x
  #     C3 <- 0.808 + 5.321*(x^0.5) - 1.231*x
  #     C4 <- -0.567 - 2.364*(x^0.5) + 0.566*x
  #   }
  #   if(x >= 2){
  #     C1 <- 1.089 + 0.924*(x^0.5) + 0.018*x
  #     C2 <- -1.504 - 2.141*(x^0.5) - 0.047*x
  #     C3 <- 2.486 + 2.289*(x^0.5) + 0.091*x
  #     C4 <- 1.056 - 1.104*(x^0.5) - 0.059*x
  #   }
  #
  #   Ktu <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  #
  #   z <- dimensions$alp
  #   CC1 <- 0.2026*(z^0.5) - 0.06620*z + 0.00281*(z^1.5)
  #   CC2 <- -0.2226*(z^0.5) - 0.07814*z + 0.002477*(z^1.5)
  #   CC3 <- 1 + 0.0298*(z^0.5) - 0.01485*z + 0.000151*(z^1.5)
  #
  #   Kt <- CC1 + CC2*(Ktu^0.5) + CC3*Ktu
  # }
  #
  # #S5. Shaft with a Shoulder Fillet, Tension (input big diameter D, radius r, little diameter d)
  # #Chart 3.4 in Peterson's Stress Concentration Factors
  # if(geometry == "shaft_fillet_tens" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
  #   t <- (dimensions$D - dimensions$d)/2
  #   x <- t/dimensions$r
  #   if(x < 0.1 || x > 20 || dimensions$D/dimensions$d < 1){
  #     stop('D > d and 0.1 < (D - d)/2r < 20 for this geometry.')
  #   }
  #
  #   if(x < 2){
  #     C1 <- 0.926 + 1.157*(x^0.5) - 0.099*x
  #     C2 <- 0.012 - 3.036*(x^0.5) + 0.961*x
  #     C3 <- -0.302 + 3.977*(x^0.5) - 1.744*x
  #     C4 <- 0.365 - 2.098*(x^0.5) + 0.878*x
  #   }
  #   if(x >= 2){
  #     C1 <- 1.2 + 0.86*(x^0.5) - 0.022*x
  #     C2 <- -1.805 - 0.346*(x^0.5) - 0.038*x
  #     C3 <- 2.198 - 0.486*(x^0.5) + 0.165*x
  #     C4 <- -0.593 - 0.028*(x^0.5) - 0.106*x
  #   }
  #
  #   y <- (2*t)/dimensions$D
  #   Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  # }
  #
  # #S6. Shaft with a Shoulder Fillet, Torsion (input big diameter D, radius r, little diameter d)
  # #Chart 3.12 in Peterson's Stress Concentration Factors
  # if(geometry == "shaft_fillet_tors" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
  #   t <- (dimensions$D - dimensions$d)/2
  #   x <- t/dimensions$r
  #   if(x < 0.25 || x > 4 || dimensions$D/dimensions$d < 1){
  #     stop('D > d and 0.25 < (D - d)/2r < 4 for this geometry.')
  #   }
  #   C1 <- 0.905 + 0.783*(x^0.5) - 0.075*x
  #   C2 <- -0.437 - 1.969*(x^0.5) + 0.553*x
  #   C3 <- 1.557 + 1.073*(x^0.5) - 0.578*x
  #   C4 <- -1.061 + 0.171*(x^0.5) + 0.086*x
  #
  #   y <- (2*t)/dimensions$D
  #   Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  # }
  # if(geometry == "shaft_regU_groove_tors" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$t) == 1){
  #   x <- dimensions$t/dimensions$r
  #   if(x < 0.25 || x > 50){
  #     stop('0.25 < t/r < 50 for this geometry.')
  #   }
  #   y <- 2*dimensions$t/dimensions$D
  #
  #   if(x == 1){
  #     C1 <- 2
  #     C2 <- -3.555
  #     C3 <- 4.898
  #     C4 <- -2.365
  #   }
  #   if(x !=1 || x < 2){
  #     C1 <- 0.966 + 1.056*(x^0.5) - 0.022*x
  #     C2 <- -0.192 - 4.037*(x^0.5) + 0.674*x
  #     C3 <- 0.808 + 5.321*(x^0.5) - 1.231*x
  #     C4 <- -0.567 - 2.364*(x^0.5) + 0.566*x
  #   }
  #   if(x >= 2){
  #     C1 <- 1.089 + 0.924*(x^0.5) + 0.018*x
  #     C2 <- -1.504 - 2.141*(x^0.5) - 0.047*x
  #     C3 <- 2.486 + 2.289*(x^0.5) + 0.091*x
  #     C4 <- 1.056 - 1.104*(x^0.5) - 0.059*x
  #   }
  #
  #   Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  # }
  # # =========================================
  # # Colin's Contributions to SCF (test now Sort later)
  # # =========================================
  # # Stepped flat tension bar with shoulder fillets, Chart 3.2a
  # if(geometry == "filled_flat_bar" && length(dimensions$L) == 1 && length(dimensions$H) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1&& length(dimensions$t) == 1){
  #   if(dimensions$L/dimensions$H >= -1.89*((dimensions$r/dimensions$d)-0.15)+5.5){
  #     stop('check dimensional requirements')
  #   }
  #   if(0.1<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=2){
  #     C1 <- 1.006 + 1.008*sqrt(dimensions$t/dimensions$r) - 0.044*(dimensions$t/dimensions$r)
  #     C2 <- -0.115 - 0.584*sqrt(dimensions$t/dimensions$r) +0.315*(dimensions$t/dimensions$r)
  #     C3 <- 0.245 - 1.006*sqrt(dimensions$t/dimensions$r) - 0.257*(dimensions$t/dimensions$r)
  #     C4 <- -0.135 + 0.582*sqrt(dimensions$t/dimensions$r) - 0.017*(dimensions$t/dimensions$r)
  #   }
  #   if(2<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=20){
  #     C1 <- 1.020 + 1.009*sqrt(dimensions$t/dimensions$r) - 0.048*(dimensions$t/dimensions$r)
  #     C2 <- -0.065 - 0.165*sqrt(dimensions$t/dimensions$r) - 0.007*(dimensions$t/dimensions$r)
  #     C3 <- -3.495 + 1.266*sqrt(dimensions$t/dimensions$r) - 0.016*(dimensions$t/dimensions$r)
  #     C4 <- 3.505 - 2.109*sqrt(dimensions$t/dimensions$r) + 0.069*(dimensions$t/dimensions$r)
  #   }
  #   x <- (2*dimensions$t)/dimensions$H
  #   Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  # }
  #
  # # Stepped tension bar of circular cross section with shoulder fillet, Chart 3.4
  # if(geometry == "filleted_circular_bar" && length(dimensions$D) == 1 && length(dimensions$d) == 1 && length(dimensions$r) == 1 && length(dimensions$t) == 1){
  #   if (loadtype==1){
  #     if(0.1<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=2){
  #       C1 <- 0.926 + 1.157*sqrt(dimensions$t/dimensions$r) - 0.099*(dimensions$t/dimensions$r)
  #       C2 <- 0.012 - 3.036*sqrt(dimensions$t/dimensions$r) + 0.961*(dimensions$t/dimensions$r)
  #       C3 <- -0.302 + 3.977*sqrt(dimensions$t/dimensions$r) - 1.744*(dimensions$t/dimensions$r)
  #       C4 <- 0.365 - 2.098*sqrt(dimensions$t/dimensions$r) + 0.878*(dimensions$t/dimensions$r)
  #     }
  #     if(2<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=20){
  #       C1 <- 1.200 + 0.860*sqrt(dimensions$t/dimensions$r) - 0.022*(dimensions$t/dimensions$r)
  #       C2 <- -1.805 - 0.346*sqrt(dimensions$t/dimensions$r) - 0.038*(dimensions$t/dimensions$r)
  #       C3 <- 2.198 - 0.486*sqrt(dimensions$t/dimensions$r) + 0.165*(dimensions$t/dimensions$r)
  #       C4 <- -0.593 - 0.028*sqrt(dimensions$t/dimensions$r) - 0.106*(dimensions$t/dimensions$r)
  #     }
  #   }
  #
  #   if (loadtype==2){
  #     if(0.1<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=2){
  #       C1 <- 0.947 + 1.206*sqrt(dimensions$t/dimensions$r) - 0.131*(dimensions$t/dimensions$r)
  #       C2 <- 0.022 - 3.405*sqrt(dimensions$t/dimensions$r) + 0.915*(dimensions$t/dimensions$r)
  #       C3 <- 0.869 + 1.777*sqrt(dimensions$t/dimensions$r) - 0.555*(dimensions$t/dimensions$r)
  #       C4 <- -0.810 + 0.422*sqrt(dimensions$t/dimensions$r) - 0.260*(dimensions$t/dimensions$r)
  #     }
  #     if(2<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=20){
  #       C1 <- 1.232 + 0.832*sqrt(dimensions$t/dimensions$r) - 0.008*(dimensions$t/dimensions$r)
  #       C2 <- -3.813 + 0.968*sqrt(dimensions$t/dimensions$r) - 0.260*(dimensions$t/dimensions$r)
  #       C3 <- 7.423 - 4.868*sqrt(dimensions$t/dimensions$r) + 0.869*(dimensions$t/dimensions$r)
  #       C4 <- -3.839 + 3.070*sqrt(dimensions$t/dimensions$r) - 0.600*(dimensions$t/dimensions$r)
  #     }
  #   }
  #   x <- (2*dimensions$t)/dimensions$D
  #   Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  # }
  #
  # #	Filleted bars in bending, Chart 3.8a
  # if(geometry == "filleted_bar" && length(dimensions$L) == 1 && length(dimensions$H) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1 && length(dimensions$t) == 1 && length(dimensions$h) == 1){
  #   if(dimensions$L/dimensions$H >= -2.05*((dimensions$r/dimensions$d)-0.025)+2){
  #     stop('check dimensional requirements')
  #   }
  #   if(0.1<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=2){
  #     C1 <- 1.006 + 0.967*sqrt(dimensions$t/dimensions$r) + 0.013*(dimensions$t/dimensions$r)
  #     C2 <- -0.270 - 2.372*sqrt(dimensions$t/dimensions$r) + 0.708*(dimensions$t/dimensions$r)
  #     C3 <- 0.662 + 1.157*sqrt(dimensions$t/dimensions$r) - 0.908*(dimensions$t/dimensions$r)
  #     C4 <- -0.405 + 0.249*sqrt(dimensions$t/dimensions$r) - 0.200*(dimensions$t/dimensions$r)
  #   }
  #   if(2<=(dimensions$t/dimensions$r) && (dimensions$t/dimensions$r)<=20){
  #     C1 <- 1.058 + 1.002*sqrt(dimensions$t/dimensions$r) - 0.038*(dimensions$t/dimensions$r)
  #     C2 <- -3.652 + 1.639*sqrt(dimensions$t/dimensions$r) - 0.436*(dimensions$t/dimensions$r)
  #     C3 <- 6.170 - 5.687*sqrt(dimensions$t/dimensions$r) + 1.175*(dimensions$t/dimensions$r)
  #     C4 <- -2.558 + 3.046*sqrt(dimensions$t/dimensions$r) - 0.701*(dimensions$t/dimensions$r)
  #   }
  #   x <- (2*dimensions$t)/dimensions$H
  #   Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  # }
  #
  # #	Tension of a thin semi-infinite element with circular hole near edge, Chart 4.2
  # if(geometry == "thin_element_edge_hole" && length(dimensions$a) == 1 && length(dimensions$c) == 1 && length(dimensions$h) == 1){
  #   x <- (dimensions$a/dimensions$c)
  #   K_tga <- 0.99619 - 0.43879*x - 0.0613028*(x^2) - 0.48941*(x^3)
  #   K_tgb <- 3.0004 + 0.083503*x + 7.3417*(x^2) - 38.046*(x^3) + 106.037*(x^4) - 130.133*(x^5) + 65.065*(x^6)
  #   K_tgc <- 2.9943 + 0.54971*x - 2.32876*(x^2) + 8.9718*(x^3) - 13.344*(x^4) + 7.1452*(x^5)
  #
  #   #this function does not calculate the K_tn, sigma and sigma_B are required for that calculation
  #   Kt <- list(K_tga, K_tgb, K_tgc)
  # }
  #
  # #	Tension of a finite-width element having an eccentrically located circular hole, Chart 4.3
  # if(geometry == "thin_element_eccentric_hole" && length(dimensions$a) == 1 && length(dimensions$c) == 1 && length(dimensions$e) == 1 && length(dimensions$h) == 1){
  #   x <- (dimensions$a/dimensions$c)
  #   y <- (dimensions$c/dimensions$e)
  #   C1 <- 2.9969 - 0.0090*y + 0.01338*(y^2)
  #   C2 <- 0.1217 + 0.5180*y - 0.5297*(y^2)
  #   C3 <- 0.5565 + 0.7215*y + 0.6153*(y^2)
  #   C4 <- 4.082 + 6.0146*y - 3.9815*(y^2)
  #   K_tg <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  #
  #   C1 <- 2.989 - 0.0064*y
  #   C2 <- -2.872 + 0.095*y
  #   C3 <- 2.348 + 0.196*y
  #   K_tn <- C1 + C2*x + C3*(x^2)
  #
  #   Kt <- list(K_tn, K_tg)
  # }
  #
  # #	Circular hole in a cylindrical shell in tension, Chart 4.4
  # if(geometry == "cylindrical_shell_circular_hole" && length(dimensions$a) == 1 && length(dimensions$R) == 1 && length(dimensions$h) == 1){
  #   B <- (((3*(1-(1/3)^2))^0.25)/2) * (dimensions$a/(sqrt(dimensions$R*dimensions$h)))
  #   x <- (dimensions$h/dimensions$R)
  #   if(x>0.02){
  #     stop('check dimensional requirements')
  #   }
  #   if(loadtype==2) {
  #     Ktg <- 2.9296 + 0.9865*B + 0.5576*(B^2) * 0.5576*(B^3) + 0.01952*(B^4)
  #     Ktn <- Ktg*(1 - (dimensions$a/(Pi*dimensions$R)))
  #   }
  #   if(loadtype==1){
  #     C1 <- 2.9127 - 3.4614*x + 277.38*(x^2)
  #     C2 <- 1.3633 - 1.9581*x + 1124.24*(x^2)
  #     C3 <- 1.3365 - 174.54*x + 21452.3*(x^2) - 683125*(x^3)
  #     C4 <- -0.5115 + 13.918*x - 335.338*(x^2)
  #     C5 <- 0.06154 - 1.707*x + 34.614*(x^2)
  #     Ktn <- C1 + C2*B + C3*(B^2) + C4*(B^3) + C4*(B^4)
  #   }
  #   Kt <- Ktn
  # }
  #
  # #	Circular hole in a cylindrical shell with internal pressure, Chart 4.5
  # if(geometry == "cylindrical_shell_pressurized" && length(dimensions$a) == 1 && length(dimensions$R) == 1 && length(dimensions$h) == 1){
  #   B <- (((3*(1-(1/3)^2))^0.25)/2) * (dimensions$a/(sqrt(dimensions$R*dimensions$h)))
  #
  #   if(loadtype==2) {
  #     if(B<=0.5){
  #       Kt <- 2.601 + 0.9147*B + 2.5*(B^2) + 30.556*(B^3) - 41.667*(B^4)
  #     }
  #     else if(B>0.5 && B<=3.14){
  #       Kt <- 1.392 + 7.394*B - 0.908*(B^2) + 0.4158*(B^3) - 0.06115*(B^4)
  #     }
  #     else {
  #       stop('Beta requirements not met')
  #     }
  #   }
  #
  #   if(loadtype==1){
  #     if(B<=2){
  #       Kt <- 2.5899 + 0.8002*B + 4.0112*(B^2) - 1.8235*(B^3) + 0.375*(B^4)
  #     }
  #     else if(B>2 && B<=4){
  #       Kt <- 8.3065 - 7.1716*B + 6.70*(B^2) - 1.35*(B^3) + 0.1056*(B^4)
  #     }
  #     else {
  #       stop('Beta requirements not met')
  #     }
  #   }
  # }
  #
  # #	Pressurized spherical shell with elliptical hole, Chart 4.6
  # if(geometry == "cylindrical_shell_elliptical_hole" && length(dimensions$a) == 1 && length(dimensions$b) == 1 && length(dimensions$R) == 1 && length(dimensions$h) == 1){
  #   x <- dimensions$a/dimensions$R
  #   y <- sqrt(dimensions$R/dimensions$h)
  #   ab <- dimensions$b/dimensions$a
  #
  #   C1 <- -1.9869 + 5.3403*ab - 1.556*(ab^2)
  #   C2 <- 5.4355 - 6.75*ab + 4.993*(ab^2)
  #   C3 <- -7.8057 + 13.2508*ab - 5.8544*(ab^2)
  #   C4 <- 1.9069 - 3.3306*ab + 1.4238*(ab^2)
  #
  #   Kt <- C1 + C2*x*y + C3*((x*y)^2) + C4*((x*y)^2)
  # }
  #
  # #	Tension of infinite panel with two circular holes, Chart 4.21
  # if(geometry == "infinite_panel_two_holes" && length(dimensions$d) == 1 && length(dimensions$l) == 1){
  #   x <- dimensions$d/dimensions$l
  #   Kt <- 3.000 - 0.712*x + 0.271*(x^2)
  # }
  #
  # #	Torsion of a splined shaft without a mating member, Chart 5.4
  # if(geometry == "splined_shaft" && length(dimensions$d) == 1 && length(dimensions$r) == 1){
  #   if(loadtype != 3){
  #     stop('the "splined_shaft" function is only available for torsional loading')
  #   }
  #   x <- (10*dimensions$r)/dimensions$d
  #   Kt <- 6.083 - 14.775*x + 18.25*(x^2)
  # }

  return(Kt)
}
