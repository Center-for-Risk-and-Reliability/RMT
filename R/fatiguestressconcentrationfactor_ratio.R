# Fatigue Stress Concentration Factor Ratio Calculator
# Developed by Reuel Smith, 2022

fatigue.stress.concentration.factor.ratio <- function(Su,material,stressunits){
  # DATA:
  Kfrat_dat <- c(0.1115977,0.1460526,0.1805075,0.2063910,0.2322368,0.2710526,0.3098496,0.3400376,0.3659774,
                 0.4004699,0.4263722,0.4566165,0.4868609,0.5257331,0.5603008,0.5949060,0.6294925,0.6555263,
                 0.6815602,0.7118233,0.7464662,0.7639098,0.7856767)
  Steel_Su_dat <- c(73.6842,81.57895,89.47368,97.36842,103.50877,114.91228,125.43860,134.21053,144.73684,
                    154.38596,163.15789,174.56140,185.96491,200.00000,213.15789,228.07018,242.10526,257.01754,
                    271.92982,284.21053,300.87719,314.91228,330.70175)
  Al_Su_dat <- c(24.31579,  26.92105,  29.52632,  32.13158,  34.19298,  38.07018,  41.64912,  44.63158,  48.21053,
                 51.49123,  54.47368,  58.35088,  62.22807,  67.00000,  71.34211,  76.26316,  80.89474,  85.81579,
                 90.73684,  94.78947, 100.28947, 104.92105, 110.13158)
  Mag_Su_dat <- c(16.21053, 17.94737, 19.68421, 21.42105, 22.77193, 25.28070, 27.59649, 29.52632, 31.84211, 33.96491,
                  35.89474, 38.40351, 40.91228, 44.00000, 46.89474, 50.17544, 53.26316, 56.54386, 59.82456, 62.52632,
                  66.19298, 69.28070, 72.75439)


  # Check units.  The coefficients were computed with Su in ksi so the Su given will
  # need to be translated accordingly.
  if(missing(stressunits) || stressunits == 1){
    Su <- Su/6.8947572932
  }
  if(stressunits == 2){
    Su <- Su
  }

  # Only computes for steels, aluminums, and magnesiums for the moment.  Will need
  # further study for additional material relations.
  if(material == "Steel"){
    params<-lm(Kfrat_dat ~ 0 + poly(Steel_Su_dat,4,raw=TRUE))
    # C1 <- -0.000194242189507
    # C2 <- 3.59976864433E-5
    # C3 <- -1.44791840544E-7
    # C4 <- 1.8045624857E-10
  }
  if(material == "Aluminum"){
    params<-lm(Kfrat_dat ~ 0 + poly(Al_Su_dat,4,raw=TRUE))
    # C1 <- -0.000111871809443
    # C2 <- 3.05988868369E-4
    # C3 <- -3.70035054616E-6
    # C4 <- 1.3855334474E-8
  }
  if(material == "Magnesium"){
    params<-lm(Kfrat_dat ~ 0 + poly(Mag_Su_dat,4,raw=TRUE))
    # C1 <- -0.000882919043214
    # C2 <- 7.43753852135E-4
    # C3 <- -1.35980316063E-5
    # C4 <- 7.7033778674E-8
  }

  C1 <- summary(params)$coefficients[1,1]
  C2 <- summary(params)$coefficients[2,1]
  C3 <- summary(params)$coefficients[3,1]
  C4 <- summary(params)$coefficients[4,1]
  Kfratio <- sum(C1*Su,C2*(Su^2),C3*(Su^3),C4*(Su^4))
  if(Kfratio > 1){
    Kfratio <- 1
  }
  return(Kfratio)
}
