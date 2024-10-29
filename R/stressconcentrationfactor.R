# Stress Concentration Factor Calculator
# Developed by Reuel Smith, 2022-2024

stress.concentration.factor <- function(dimensions,geometry,loadtype = 1){
  # Calculates stress concentration factor Kt based on Peterson's Stress Concentration Factors
  # Default for 'loadtype' will be tension (1).  Can also have bending (2) and traverse bending (3)
  # =========================================
  # Rectangular Bar or Plate options
  # =========================================
  # RB1. Rectangular Bar with a Semi-circle Edge Notch (input width W and radius r)
  # Chart 2.9, 2.30a in Peterson's Stress Concentration Factors
  if(geometry == "rect_1semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    # if(dimensions$r/dimensions$W >= 0.3){
    #   stop('r/W must be less than 0.3 for this geometry.')
    # }
    x <- dimensions$r/dimensions$W

    if(loadtype == 1){
      # Tension
      Kt <- 3.065 - 8.871*x + 14.036*(x^2) - 7.219*(x^3)
    }
    if(loadtype == 2){
      # Bending
      Kt <- 3.065 - 6.643*x + 0.205*(x^2) - 4.004*(x^3)
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }
  }
  # RB2. Rectangular Bar with Opposite Semi-circle Edge Notches (input width W and radius r)
  # Chart 2.3, 2.24 in Peterson's Stress Concentration Factors
  if(geometry == "rect_2semicirc_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1){
    if(2*dimensions$r >= dimensions$W){
      stop('Radius exceeds the width of this geometry.')
    }
    # if(dimensions$r > 0.25*dimensions$W){
    #   stop('r/W must be less than 0.5 for this geometry.')
    # }

    if(loadtype == 1){
      # Tension
      x <- (2*dimensions$r)/dimensions$W
      Kt <- 3.065 - 3.472*x + 1.009*(x^2) + 0.405*(x^3)
    }
    if(loadtype == 2){
      # Bending
      x <- (2*dimensions$r)/dimensions$W
      Kt <- 3.065 - 6.637*x + 8.229*(x^2) - 3.636*(x^3)
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }
  }
  # RB3. Rectangular Bar with a U-Notch (input width W, radius r, and U depth d)
  # Chart 2.9, 2.30a in Peterson's Stress Concentration Factors
  if(geometry == "rect_1U_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    x <- dimensions$d/dimensions$W
    y <- dimensions$d/dimensions$r
    if(dimensions$d/dimensions$W >= 0.3 || dimensions$d/dimensions$r < 0.5 || dimensions$d/dimensions$r > 20) {
      stop(' 0.5 < d/r < 20 for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      if(dimensions$d/dimensions$r == 1){
        C1 <- 3.065
        C2 <- -8.871
        C3 <- 14.036
        C4 <- -7.219
      }
      if(y != 1 && y < 2){
        C1 <- 0.907 + 2.125*sqrt(y) + 0.023*y
        C2 <- 0.701 - 11.289*sqrt(y) + 1.708*y
        C3 <- -0.672 + 18.754*sqrt(y) - 4.046*y
        C4 <- 0.175 - 9.759*sqrt(y) + 2.365*y
      }
      if(y != 1 && y >= 2){
        C1 <- 0.953 + 2.136*sqrt(y) - 0.005*y
        C2 <- -3.255 - 6.281*sqrt(y) + 0.068*y
        C3 <- 8.203 + 6.893*sqrt(y) + 0.064*y
        C4 <- -4.851 - 2.793*sqrt(y) - 0.128*y
      }
    }
    if(loadtype == 2){
      # Bending
      if(dimensions$d/dimensions$r == 1){
        C1 <- 3.065
        C2 <- -6.643
        C3 <- 0.205
        C4 <- -4.004
      }
      if(y != 1 && y < 2){
        C1 <- 1.795 + 1.481*y + 0.211*(y^2)
        C2 <- -3.544 - 3.677*y + 0.578*(y^2)
        C3 <- 5.459 + 3.691*y - 0.565*(y^2)
        C4 <- -2.678 - 1.531*y + 0.205*(y^2)
      }
      if(y != 1 && y >= 2){
        C1 <- 2.966 + 0.502*y - 0.009*(y^2)
        C2 <- -6.475 - 1.126*y + 0.019*(y^2)
        C3 <- 8.023 + 1.253*y - 0.020*(y^2)
        C4 <- -3.572 - 0.634*y + 0.010*(y^2)
      }
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }
    Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  }
  # RB4. Rectangular Bar with Opposite Edge U-Notches (input width W, radius r, and U depth d)
  # Chart 2.4, 2.25 in Peterson's Stress Concentration Factors
  if(geometry == "rect_2U_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    # if(dimensions$d/dimensions$W > 0.25 || dimensions$d/dimensions$W < 0.02 || dimensions$d/dimensions$r < 0.1 || dimensions$d/dimensions$r > 50) {
    #   stop('0.02 < d/W < 0.25 and 0.1 < d/r < 50 for this geometry.')
    # }
    if(dimensions$d/dimensions$r < 0.1 || dimensions$d/dimensions$r > 50) {
      stop('0.1 < d/r < 50 for this geometry.')
    }
    x <- (2*dimensions$d)/dimensions$W
    y <- dimensions$d/dimensions$r

    if(loadtype == 1){
      # Tension
      # if(dimensions$d/dimensions$W > 0.25 || dimensions$d/dimensions$W < 0.02 || y < 0.1 || y > 50) {
      #   stop('0.02 < d/W < 0.25 and 0.1 < d/r < 50 for this geometryy.')
      # }
      if(y == 1){
        C1 <- 3.065
        C2 <- -3.472
        C3 <- 1.009
        C4 <- 0.405
      }
      if(y != 1 && y < 2){
        C1 <- 0.955 + 2.169*sqrt(y) + 0.081*y
        C2 <- -1.557 - 4.046*sqrt(y) + 1.032*y
        C3 <- 4.013 + 0.424*sqrt(y) - 0.748*y
        C4 <- -2.461 + 1.538*sqrt(y) - 0.236*y
      }
      if(y != 1 && y >= 2){
        C1 <- 1.037 + 1.991*sqrt(y) + 0.002*y
        C2 <- -1.886 - 2.181*sqrt(y) - 0.048*y
        C3 <- 0.649 + 1.086*sqrt(y) + 0.142*y
        C4 <- 1.218 - 0.922*sqrt(y) - 0.086*y
      }
    }
    if(loadtype == 2){
      # Bending
      if(dimensions$d/dimensions$W > 0.25 || dimensions$d/dimensions$W < 0.02 || y < 0.1 || y > 50) {
        stop('0.02 < d/W < 0.25 and 0.1 < d/r < 50 for this geometry.')
      }
      if(y == 1){
        C1 <- 3.065
        C2 <- -6.637
        C3 <- 8.229
        C4 <- -3.636
      }
      if(y != 1 && y < 2){
        C1 <- 1.024 + 2.092*sqrt(y) - 0.051*y
        C2 <- -0.630 - 7.194*sqrt(y) + 1.288*y
        C3 <- 2.117 + 8.574*sqrt(y) - 2.160*y
        C4 <- -1.420 - 3.494*sqrt(y) + 0.932*y
      }
      if(y != 1 && y >= 2){
        C1 <- 1.113 + 1.957*sqrt(y)
        C2 <- -2.579 - 4.017*sqrt(y) - 0.013*y
        C3 <- 4.100 + 3.922*sqrt(y) + 0.083*y
        C4 <- -1.528 - 1.893*sqrt(y) - 0.066*y
      }
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }
    Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  }
  # RB5. Rectangular Bar with a V-Notch (input width W, radius r, V depth d, and V angle alp)
  # Chart 2.7? in Peterson's Stress Concentration Factors

  # RB6. Rectangular Bar with Opposite Edge V-Notches (input width W, radius r, V-depth d, and V-angle alp)
  # Chart 2.7 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "rect_2V_edge" && length(dimensions$W) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1 && length(dimensions$alp) == 1){
    if(dimensions$alp > 150  || !(setequal((2*dimensions$d/dimensions$w),0.4) || !(setequal((2*dimensions$d/dimensions$w),2/3))) || dimensions$d/dimensions$r <= 1 || dimensions$d/dimensions$r > 50){
      stop('alpha < 150, 2d/w = 0.4 or 0.667, and 1 < d/r < 50 for this geometry.')
    }

    if(loadtype == 1){
      # Tension
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
      Ktu <- C1 + C2*x + C3*(x^2) + C4*(x^3)

      if(2*dimensions$d/dimensions$w == 0.4){
        if(Ktu > 3.5 || Ktu < 1.6){
          stop('1.6 < Ktu < 3.5 for these parameters.')
        }
        if(dimensions$alp < 90){
          Kt <- Ktu
        }

        CC1 <- 5.294 - 0.1225*dimensions$alp + 0.000523*(dimensions$alp^2)
        CC2 <--5.0002 + 0.1171*dimensions$alp - 0.000434*(dimensions$alp^2)
        CC3 <-1.423 - 0.00197*dimensions$alp - 0.000004*(dimensions$alp^2)

        Kt <- CC1 + CC2*(Ktu^0.5) + CC3*Ktu
      }

      if(2*dimensions$d/dimensions$w == 2/3){
        if(Ktu > 2.8 || Ktu < 1.6){
          stop('1.6 < Ktu < 2.8 for these parameters.')
        }
        if(dimensions$alp < 60){
          Kt <- Ktu
        }
        CC1 <- -10.01 + 0.1534*dimensions$alp - 0.000647*(dimensions$alp^2)
        CC2 <- 13.60 - 0.2140*dimensions$alp + 0.000973*(dimensions$alp^2)
        CC3 <- -3.781 + 0.07875*dimensions$alp - 0.000392*(dimensions$alp^2)

        Kt <- CC1 + CC2*(Ktu^0.5) + CC3*Ktu
      }
    }
    if(loadtype == 2){
      # Bending
      stop('Bending calculation unavailable for this geometry')
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }
  }

  # RB7. Plate with a Circular Symmetric Hole (input width W and diameter d)
  # Chart 4.1 in Peterson's Stress Concentration Factors
  if(geometry == "rect_circ_symm_hole" && length(dimensions$W) == 1 && length(dimensions$d) == 1){
    if(dimensions$d >= 0.9*dimensions$W){
      stop('Measurement for d must be less than 0.9*W for this geometry.')
    }
    x <- 1 - dimensions$d/dimensions$W
    Kt <- 2 + 0.284*x - 0.6*(x^2) + 1.32*(x^3)
  }
  # RB8. Plate with a Elliptical Symmetric Hole (input width W, minor diameter a, and major diameter d)
  # Chart 4.51 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "rect_ellips_symm_hole" && length(dimensions$W) == 1 && length(dimensions$a) == 1 && length(dimensions$d) == 1){
    if(dimensions$d/dimensions$a < 1 || dimensions$d/dimensions$a > 8){
      stop('d/a must be between 1 and 8 for this geometry.')
    }
    C1 <- 1.109 - 0.188*sqrt(dimensions$d/dimensions$a) + 2.086*(dimensions$d/dimensions$a)
    C2 <- -0.486 + 0.213*sqrt(dimensions$d/dimensions$a) - 2.588*(dimensions$d/dimensions$a)
    C3 <- 3.816 - 5.510*sqrt(dimensions$d/dimensions$a) + 4.638*(dimensions$d/dimensions$a)
    C4 <- -2.438 + 5.485*sqrt(dimensions$d/dimensions$a) - 4.126*(dimensions$d/dimensions$a)

    x <- dimensions$d/dimensions$W
    Kt <- C1 + C2*x + C3*(x^2) + C4*(x^3)
  }

  # RB9. Plate with a Slot Hole (input width W, slot distance d, and radius r)
  # Chart 4.1 in Peterson's Stress Concentration Factors

  # RB10. Plate with a Circular Offset Hole (input width W, edge to hole center distance a, and radius r)
  # Chart 4.2 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "plate_circ_offset_hole" && length(dimensions$w) == 1 && length(dimensions$a) == 1 && length(dimensions$r) == 1){
    if(dimensions$r/dimensions$a > 0.8 || (dimensions$w-dimensions$a)/dimensions$a < 1){
      stop('r/a must be less than 0.8 for this geometry and (w-a)/a < 1.')
    }
    x <- (dimensions$a/(dimensions$w-dimensions$a))
    C1 <- 2.989 - 0.0064*x
    C2 <- -2.872 + 0.095*x
    C3 <- 2.348 + 0.196*x

    y <- (dimensions$r/dimensions$a)
    Kt <- C1 + C2*y + C3*(y^2)
  }

  # RB11. Plate with a Deep Hyperbolic Notch (input edge to hyperbola distance d and thickness t)
  # Chart 4.2 in Peterson's Stress Concentration Factors

  # RB12. Plate with Two Deep Hyperbolic Notches (input edge to hyperbola distance d and thickness t)
  # Chart 2.1 in Peterson's Stress Concentration Factors

  # RB13. Rectangular Bar with Infinite Rows of Semicircular Notches (input width w, spacing b, radius r)
  # Chart 2.12 in Peterson's Stress Concentration Factors  (CAMILLE)
  if(geometry == "rect_inf_semicirc_edge" && length(dimensions$w) == 1 && length(dimensions$b) == 1 && length(dimensions$r) == 1){
    if(2*dimensions$r/dimensions$w > 0.4 || 2*dimensions$r/dimensions$b > 1){
      stop('2r/w must be less than 0.4 for this geometry and 2r/b must be < 1.')
    }
    x <- ((2*dimensions$r)/dimensions$w)
    C1 <- 3.1055 - 3.2487*x + 0.8522*(x^2)
    C2 <- - 1.437 + 10.5053*x - 8.7547*(x^2) + 19.6273*(x^3)
    C3 <- -1.6753 - 14.0851*x + 43.675*(x^2)
    C4 <- 1.7207 + 5.7974*x - 27.7453*(x^2) + 6.0444*(x^3)

    y <- ((2*dimensions$r)/dimensions$b)
    Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  }

  # RB15. Rectangular Bar with Double Fillets (input width W, fillet spacing L, inner width d, radius r)
  # Chart 3.2a and 3.8a in Peterson's Stress Concentration Factors  (COLIN)
  if(geometry == "rect_double_fillet" && length(dimensions$W) == 1 && length(dimensions$L) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    t <- (dimensions$W-dimensions$d)/2
    x <- t/dimensions$r
    if(loadtype == 1 && dimensions$L/dimensions$W >= -1.89*((dimensions$r/dimensions$d)-0.15)+5.5){
      stop('L/W > -1.89(r/d - 0.15) + 5.5 for this geometry')
    }
    if(loadtype == 2 && dimensions$L/dimensions$W >= -2.05*((dimensions$r/dimensions$d)-0.025)+2){
      stop('L/W > -2.05(r/d - 0.025) + 2.0 for this geometry')
    }
    if(x < 0.1 || x > 20){
      stop('0.1 < a/r < 20 for this geometry')
    }
    if(loadtype == 1){
      # Tension
      if(0.1<=x && x<=2){
        C1 <- 1.006 + 1.008*sqrt(x) - 0.044*x
        C2 <- -0.115 - 0.584*sqrt(x) +0.315*x
        C3 <- 0.245 - 1.006*sqrt(x) - 0.257*x
        C4 <- -0.135 + 0.582*sqrt(x) - 0.017*x
      }
      if(x>=2 && x<=20){
        C1 <- 1.020 + 1.009*sqrt(x) - 0.048*x
        C2 <- -0.065 - 0.165*sqrt(x) - 0.007*x
        C3 <- -3.495 + 1.266*sqrt(x) - 0.016*x
        C4 <- 3.505 - 2.109*sqrt(x) + 0.069*x
      }
    }
    if(loadtype == 2){
      # Bending
      if(0.1<=x && x<=2){
        C1 <- 1.006 + 0.967*sqrt(x) + 0.013*x
        C2 <- -0.270 - 2.372*sqrt(x) + 0.708*x
        C3 <- 0.662 + 1.157*sqrt(x) - 0.908*x
        C4 <- -0.405 + 0.249*sqrt(x) - 0.200*x
      }
      if(x>=2 && x<=20){
        C1 <- 1.058 + 1.002*sqrt(x) - 0.038*x
        C2 <- -3.652 + 1.639*sqrt(x) - 0.436*x
        C3 <- 6.170 - 5.687*sqrt(x) + 1.175*x
        C4 <- -2.558 + 3.046*sqrt(x) - 0.701*x
      }
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }

    y <- (2*t)/dimensions$W
    Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  }

  # =========================================
  # Round or Shaft Options
  # =========================================
  #RS1. Round Shaft with a Single Fillet (input outer diameter D, inner diameter d, radius r)
  #Chart 3.4, 3.10, and 3.12 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "shaft_fillet_single" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    t <- (dimensions$D - dimensions$d)/2
    x <- t/dimensions$r
    if((x < 0.1 || x > 20 || dimensions$D/dimensions$d < 1) && (loadtype == 1 || loadtype == 2)){
      stop('D > d and 0.1 < (D - d)/2r < 20 for this geometry.')
    }
    if((x < 0.25 || x > 4 || dimensions$D/dimensions$d < 1) && loadtype == 3){
      stop('D > d and 0.25 < (D - d)/2r < 4 for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      if(x < 2){
        C1 <- 0.926 + 1.157*(x^0.5) - 0.099*x
        C2 <- 0.012 - 3.036*(x^0.5) + 0.961*x
        C3 <- -0.302 + 3.977*(x^0.5) - 1.744*x
        C4 <- 0.365 - 2.098*(x^0.5) + 0.878*x
      }
      if(x >= 2){
        C1 <- 1.2 + 0.86*(x^0.5) - 0.022*x
        C2 <- -1.805 - 0.346*(x^0.5) - 0.038*x
        C3 <- 2.198 - 0.486*(x^0.5) + 0.165*x
        C4 <- -0.593 - 0.028*(x^0.5) - 0.106*x
      }
    }
    if(loadtype == 2){
      # Bending
      if(x < 2){
        C1 <- 0.947 + 1.206*(x^0.5) - 0.131*x
        C2 <- 0.022 - 3.405*(x^0.5) + 0.915*x
        C3 <- 0.869 + 1.777*(x^0.5) - 0.555*x
        C4 <- -0.810 + 0.422*(x^0.5) - 0.260*x
      }
      if(x >= 2){
        C1 <- 1.232 + 0.832*(x^0.5) - 0.008*x
        C2 <- -3.813 - 0.968*(x^0.5) - 0.260*x
        C3 <- 7.423 - 4.868*(x^0.5) + 0.869*x
        C4 <- -3.839 + 3.070*(x^0.5) - 0.600*x
      }
    }
    if(loadtype == 3){
      # Torsion
      C1 <- 0.905 + 0.783*(x^0.5) - 0.075*x
      C2 <- -0.437 - 1.969*(x^0.5) + 0.553*x
      C3 <- 1.557 + 1.073*(x^0.5) - 0.578*x
      C4 <- -1.061 + 0.171*(x^0.5) + 0.086*x
    }

    y <- (2*t)/dimensions$D
    Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
  }

  #RS4. Round Shaft with a U-Shaped or Semicircular Groove (input diameter D, radius of groove r, depth of groove t)
  #Chart 2.47 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "shaft_reg_U_groove" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$t) == 1){
    x <- dimensions$t/dimensions$r
    if(x < 0.25 || x > 50){
      stop('0.25 < t/r < 50 for this geometry.')
    }
    if(loadtype == 1 && dimensions$t == dimensions$r){
      # Tension
      C1 <- 3.02355888095581
      C2 <- -4.79623997372538
      C3 <- 4.64010298684958
      C4 <- -1.00930131841631
    }
    if(loadtype == 2){
      # Bending
      if(x == 1){
        C1 <- 3.032
        C2 <- -7.431
        C3 <- 10.390
        C4 <- -5.009
      }
      if(x !=1 && (x > 0.25|| x < 2)){
        C1 <- 0.594 + 2.958*(x^0.5) - 0.520*x
        C2 <- 0.422 - 10.545*(x^0.5) + 2.692*x
        C3 <- 0.501 + 14.375*(x^0.5) - 4.486*x
        C4 <- -0.613 - 6.573*(x^0.5) + 2.177*x
      }
      if(x >= 2){
        C1 <- 0.965 + 1.926*(x^0.5)
        C2 <- -2.773 - 4.414*(x^0.5) - 0.017*x
        C3 <- 4.785 + 4.681*(x^0.5) + 0.096*x
        C4 <- -1.995 - 2.241*(x^0.5) - 0.074*x
      }
    }
    if(loadtype == 3){
      # Torsion
      if(x == 1){
        C1 <- 2
        C2 <- -3.555
        C3 <- 4.898
        C4 <- -2.365
      }
      if(x !=1 && (x > 0.25|| x < 2)){
        C1 <- 0.966 + 1.056*(x^0.5) - 0.022*x
        C2 <- -0.192 - 4.037*(x^0.5) + 0.674*x
        C3 <- 0.808 + 5.321*(x^0.5) - 1.231*x
        C4 <- -0.567 - 2.364*(x^0.5) + 0.566*x
      }
      if(x >= 2){
        C1 <- 1.089 + 0.924*(x^0.5) + 0.018*x
        C2 <- -1.504 - 2.141*(x^0.5) - 0.047*x
        C3 <- 2.486 + 2.289*(x^0.5) + 0.091*x
        C4 <- 1.056 - 1.104*(x^0.5) - 0.059*x
      }
    }
    if((dimensions$t == dimensions$r && loadtype == 1) || loadtype == 2 || loadtype == 3){
      y <- 2*dimensions$t/dimensions$D
      Kt <- C1 + C2*y + C3*(y^2) + C4*(y^3)
    }
    if(dimensions$t != dimensions$r && loadtype == 1){
      y <- (dimensions$D - 2*dimensions$t)/dimensions$D
      Kp <- sqrt(x)*sqrt(1/(1-y)) - 1
      Kq <- dimensions$t*sqrt(x)
      Kt <- 1 + 1/sqrt((1/(1.197*Kp))^2 + (1/(1.871*Kq))^2)
    }
  }

  #RS5. Shaft with a V-Shaped Groove, Torsion (input diameter D, radius r, V depth d, and V angle alp)
  #Chart 2.51 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "shaft_reg_V_groove" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1 && length(dimensions$alp) == 1){
    dia <- dimensions$D - (2*dimensions$d)

    if(dimensions$alp > 125 || (dimensions$alp > 90 && dimensions$r/dia > 0.01)){
      stop('alpha < 125 for this geometry, and for alpha > 90, r/(D-2d) < 0.01.')
    }

    x <- dimensions$d/dimensions$r
    if(x < 0.25 || x > 50){
      stop('0.25 < d/r < 50 for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      stop('Tension calculation unavailable for this geometry')
    }
    if(loadtype == 2){
      # Bending
      if(x == 1){
        C1 <- 3.032
        C2 <- -7.431
        C3 <- 10.390
        C4 <- -5.009
      }
      if(x !=1 && (x > 0.25|| x < 2)){
        C1 <- 0.594 + 2.958*(x^0.5) - 0.520*x
        C2 <- 0.422 - 10.545*(x^0.5) + 2.692*x
        C3 <- 0.501 + 14.375*(x^0.5) - 4.486*x
        C4 <- -0.613 - 6.573*(x^0.5) + 2.177*x
      }
      if(x >= 2){
        C1 <- 0.965 + 1.926*(x^0.5)
        C2 <- -2.773 - 4.414*(x^0.5) - 0.017*x
        C3 <- 4.785 + 4.681*(x^0.5) + 0.096*x
        C4 <- -1.995 - 2.241*(x^0.5) - 0.074*x
      }
    }
    if(loadtype == 3){
      # Torsion
      if(x == 1){
        C1 <- 2
        C2 <- -3.555
        C3 <- 4.898
        C4 <- -2.365
      }
      if(x !=1 && (x > 0.25|| x < 2)){
        C1 <- 0.966 + 1.056*(x^0.5) - 0.022*x
        C2 <- -0.192 - 4.037*(x^0.5) + 0.674*x
        C3 <- 0.808 + 5.321*(x^0.5) - 1.231*x
        C4 <- -0.567 - 2.364*(x^0.5) + 0.566*x
      }
      if(x >= 2){
        C1 <- 1.089 + 0.924*(x^0.5) + 0.018*x
        C2 <- -1.504 - 2.141*(x^0.5) - 0.047*x
        C3 <- 2.486 + 2.289*(x^0.5) + 0.091*x
        C4 <- 1.056 - 1.104*(x^0.5) - 0.059*x
      }
    }

    y <- 2*dimensions$d/dimensions$D
    Ktu <- C1 + C2*y + C3*(y^2) + C4*(y^3)

    z <- dimensions$alp
    CC1 <- 0.2026*(z^0.5) - 0.06620*z + 0.00281*(z^1.5)
    CC2 <- -0.2226*(z^0.5) - 0.07814*z + 0.002477*(z^1.5)
    CC3 <- 1 + 0.0298*(z^0.5) - 0.01485*z + 0.000151*(z^1.5)

    Kt <- CC1 + CC2*(Ktu^0.5) + CC3*Ktu
  }

  #RS8. Cylindrical Tube with a Circular Hole (input outer radius R, inner radius r, and hole diameter d)
  # Chart 4.4, 4.107, and 4.108 in Peterson's Stress Concentration Factors (COLIN)
  if(geometry == "shaft_cyl_tube_circ_hole" && length(dimensions$R) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    h <- dimensions$R - dimensions$r
    a <- 0.5*dimensions$d
    Rp <- 0.5*(dimensions$R + dimensions$r)
    B <- (((3*(1 - ((1/3)^2)))^0.25)/2)*(dimensions$a/(sqrt(Rp*h)))
    x <- h/Rp

    if(x > 0.02){
      stop('2(R - r)/(R + r) < 0.02 for this geometry.')
    }

    if(loadtype == 1){
      # Tension
      C1 <- 2.9127 - 3.4614*x + 277.38*(x^2)
      C2 <- 1.3633 - 1.9581*x + 1124.24*(x^2)
      C3 <- 1.3365 - 174.54*x + 21452.3*(x^2) - 683125*(x^3)
      C4 <- -0.5115 + 13.918*x - 335.338*(x^2)
      C5 <- 0.06154 - 1.707*x + 34.614*(x^2)
      Ktn <- C1 + C2*B + C3*(B^2) + C4*(B^3) + C4*(B^4)
    }
    if(loadtype == 2){
      # Bending
      Ktg <- 2.9296 + 0.9865*B + 0.5576*(B^2) * 0.5576*(B^3) + 0.01952*(B^4)
      Ktn <- Ktg*(1 - (dimensions$a/(pi*dimensions$R)))
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry.  Check next update.')
    }
    Kt <- Ktn
  }

  #RS9. Round Bar with as Shallow U-Shaped (dog-bone) Groove (input outer diameter D, inner diameter d, and radius of groove r)
  # Chart 2.21, 2.43, and 2.50 in Peterson's Stress Concentration Factors (CAMILLE)
  if(geometry == "shaft_wideU_groove" && length(dimensions$D) == 1 && length(dimensions$r) == 1 && length(dimensions$d) == 1){
    if(dimensions$D/dimensions$d > 1.1 || dimensions$D/dimensions$d < 1.005 || dimensions$r/dimensions$d > 1 || dimensions$r/dimensions$d < 0.3){
      stop('1.005 < D/d < 1.1 and 0.3 < r/d < 1 for this geometry.')
    }
    x <- dimensions$D/dimensions$d
    if(loadtype == 1){
      # Tension
      C1 <- -81.39 + 153.1*x - 70.49*(x^2)
      C2 <- 119.64 - 221.81*x + 101.93*(x^2)
      C3 <- -57.88 + 107.33*x - 49.34*(x^2)
    }
    if(loadtype == 2){
      # Bending
      C1 <- -39.58 + 73.22*x - 32.46*(x^2)
      C2 <- -9.477 + 29.41*x - 20.13*(x^2)
      C3 <- 82.46 - 166.96*x + 84.58*(x^2)
    }
    if(loadtype == 3){
      # Torsion
      C1 <- -35.16 + 67.57*x - 31.28*(x^2)
      C2 <- 79.13 - 148.37*x + 69.09*(x^2)
      C3 <- -50.34 + 94.67*x - 44.26*(x^2)
    }


    y <- dimensions$r/dimensions$d
    Kt <- C1 + C2*y + C3*(y^2)
  }

  # =========================================
  # Infinite/Semi-Infinite Plate Options
  # =========================================
  #	ISI2. Infinite Plate with Two Circular Holes (input distance between two holes a and radius r),
  # Chart 4.21b and 4.22 in Peterson's Stress Concentration Factors (COLIN)
  if(geometry == "inf_plate_two_circ_holes" && length(dimensions$r) == 1 && length(dimensions$a) == 1){
    x <- (2*dimensions$r)/dimensions$a
    if(x > 10){
      stop('0 < 2r/a <= 10 for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      Kt <- 3.000 - 0.712*x + 0.271*(x^2)
    }
    if(loadtype == 2){
      # Bending
      Kt <- 3.000 - 3.0018*x + 1.0099*(x^2)
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry.  Check next update.')
    }
  }

  # ISI4. Semi-Infinite Plate with a Circular Offset Hole (input edge to hole center distance a, and radius r)
  # Chart 4.2 in Peterson's Stress Concentration Factors (COLIN)
  if(geometry == "semi_inf_plate_circ_offset_hole" && length(dimensions$a) == 1 && length(dimensions$r) == 1){
    x <- (dimensions$r/dimensions$a)
    if(x >= 1){
      stop('r < a for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      K_tgb <- 3.0004 + 0.083503*x + 7.3417*(x^2) - 38.046*(x^3) + 106.037*(x^4) - 130.133*(x^5) + 65.065*(x^6)
    }

    Kt <- K_tgb*((1 - x)/sqrt(1 - x^2))
  }
  # ISI5. Infinite Plate with Elliptical hole (input major radius (a) and minor radius (b) of elliptical hole)
  # Chart 4.50 in Peterson's Stress Concentration Factors (MATT THOMAS) RCS - 09192024
  if(geometry == "inf_plate_ellips_hole" && length(dimensions$a) == 1 && length(dimensions$b) == 1){
    x <- (dimensions$a/dimensions$b)
    if(x < 0 || x > 10){
      stop('a/b must be between 0 and 10 for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      K_tg <- 1 + 2*x
    }
    if(loadtype == 2){
      # Bending
      stop('Bending calculation unavailable for this geometry')
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry.  Check next update.')
    }
    Kt <- K_tg
  }

  # ISI6. Semi-Infinite Plate with a U-shaped notch on one edge (input radius r and U-depth d)
  # Chart 2.2 in Peterson's Stress Concentration Factors (can also be for elliptical edge notches)
  if(geometry == "semi_inf_plate_1U_edge" && length(dimensions$d) == 1 && length(dimensions$r) == 1){
    x <- (dimensions$d/dimensions$r)
    if(x < 1){
      stop('r <= d for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      K_t <- 0.855 + 2.21*sqrt(x)
    }
    if(loadtype == 2){
      # Bending
      stop('Bending calculation unavailable for this geometry')
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry.')
    }
  }

  # ISI7. Semi-Infinite Plate with a elliptical-shaped notch on one edge (input smallest radius r and elliptical depth d)
  # Chart 2.2, 2.37 in Peterson's Stress Concentration Factors (can also be for elliptical edge notches)
  if(geometry == "semi_inf_plate_ellips_edge" && length(dimensions$d) == 1 && length(dimensions$r) == 1){
    x <- (dimensions$d/dimensions$r)
    if(x < 1){
      stop('r <= d for this geometry.')
    }
    if(loadtype == 1){
      # Tension
      K_t <- 0.855 + 2.21*sqrt(x)
    }
    if(loadtype == 2){
      # Bending
      K_t <- 0.998 + 0.79*sqrt(x)
    }
    if(loadtype == 3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry.')
    }
  }

  # =========================================
  # Pressurized Round or Shaft Options
  # =========================================
  #PRS1. Pressurized (Internal) Cylindrical Tube with a Circular Hole (input outer radius R, inner radius r, and hole diameter d)
  # Chart 4.5 in Peterson's Stress Concentration Factors (COLIN)
  if(geometry == "shaft_cyl_tube_circ_hole_pressurized" && length(dimensions$r) == 1 && length(dimensions$R) == 1 && length(dimensions$d) == 1){
    a <- dimensions$d/2
    Rp <- dimensions$r
    h <- dimensions$R - dimensions$r
    B <- (((3*(1-(1/3)^2))^0.25)/2) * (a/(sqrt(Rp*h)))

    if(loadtype==1){
      if(B<=2){
        Kt <- 2.5899 + 0.8002*B + 4.0112*(B^2) - 1.8235*(B^3) + 0.375*(B^4)
      }
      else if(B>2 && B<=4){
        Kt <- 8.3065 - 7.1716*B + 6.70*(B^2) - 1.35*(B^3) + 0.1056*(B^4)
      }
      else {
        stop('d/sqrt(r(R - r)) needs to be between 0 and 12.52068')
      }
    }

    if(loadtype==2) {
      if(B<=0.5){
        Kt <- 2.601 + 0.9147*B + 2.5*(B^2) + 30.556*(B^3) - 41.667*(B^4)
      }
      else if(B>0.5 && B<=3.14){
        Kt <- 1.392 + 7.394*B - 0.908*(B^2) + 0.4158*(B^3) - 0.06115*(B^4)
      }
      else {
        stop('d/sqrt(r(R - r)) needs to be between 0 and 9.828731')
      }
    }

    if(loadtype==3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }
  }

  #PRS2. Pressurized Spherical Shell with Elliptical Hole (input outer radius of sphere R, inner radius of sphere r, semi-major radius of ellipse a, and semi-minor radius of ellipse b)
  #	Chart 4.6 in Peterson's Stress Concentration Factors (COLIN)
  if(geometry == "shell_spherical_ellips_hole_pressurized" && length(dimensions$a) == 1 && length(dimensions$b) == 1 && length(dimensions$R) == 1 && length(dimensions$r) == 1){
    h <- dimensions$R - dimensions$r
    x <- dimensions$a/dimensions$r
    y <- sqrt(dimensions$r/h)
    ab <- dimensions$b/dimensions$a

    if(loadtype==1){
      # Tension
      C1 <- -1.9869 + 5.3403*ab - 1.556*(ab^2)
      C2 <- 5.4355 - 6.75*ab + 4.993*(ab^2)
      C3 <- -7.8057 + 13.2508*ab - 5.8544*(ab^2)
      C4 <- 1.9069 - 3.3306*ab + 1.4238*(ab^2)
    }

    if(loadtype==2) {
      # Bending
      stop('Bending calculation unavailable for this geometry')
    }

    if(loadtype==3){
      # Torsion
      stop('Torsion calculation unavailable for this geometry')
    }



    Kt <- C1 + C2*x*y + C3*((x*y)^2) + C4*((x*y)^2)
  }

  return(Kt)
}
