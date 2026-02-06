# Summary of the Reliability Modeling Toolkit (RMT) Workshop held in Mirimar Beach
# FL for the RAMS 2026 Conference
# ==============================================================================
# INSTALLATION INSTRUCTIONS
# ==============================================================================
# 1. Download and install latest version of R (https://www.r-project.org/ ).

# 2. RStudio recommended as an additional download for interface reasons
# (https://posit.co/download/rstudio-desktop/).  If installing on a Mac,
# use Anaconda Navigator to install RStudio.

# 3. Setup R and RStudio then install the “devtools” library
install.packages("devtools")

# 4. Install “cmdstanr” library
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# 5. Activate the “devtools” library to connect with GitHub.
# Then Install “RMT”
# If Rtools is not installed (usually the case), type the following to build from source,
devtools::install_github("Center-for-Risk-and-Reliability/RMT", build = FALSE, INSTALL_opts = "--install-tests")
# If Rtools is installed, type the following to build from source
devtools::install_github("Center-for-Risk-and-Reliability/RMT", INSTALL_opts = "--install-tests")

# 6. Finally, load RMT library in R
library(reliabilityRMT)
# Installation includes additional libraries that may need updating

# ==============================================================================
# PART 1: RMT Basics – Probability Plotting, MLE, and Bayesian Estimation
# ==============================================================================
# A. Probability Plotting
# Example 1: Data entry and probability plotting
data_RAMS.Pt1.Ex1 <- cbind(c(2350, 2560, 1980,
                             220,	250, 330, 370, 380, 460, 460, 510, 610),
                           c(rep(1,3),rep(1,9)),
                           c(rep(150+273,3),rep(200+273,9)))
probplot.logn(data_RAMS.Pt1.Ex1,"Blom","Time to Failure (hours)",stressunit1 = "K")
# Try this with the Weibull probablity plotting tool probplot.wbl

# Example 2: Data entry with right-censored data and probability plotting
# DIY Example with Right-Censored Data
data_RAMS.Pt1.Ex2 <- cbind(c(250, 460, 530, 730, 820, 970, 1530,
                             970,
                             160, 180, 290, 320, 390, 460,
                             500, 500,
                             90, 100, 150, 180, 220, 230, 230,
                             250),
                           c(rep(1,7),0,rep(1,6),rep(0,2),rep(1,7),0),
                           c(rep(200,8),rep(300,8),rep(475,8)))
probplot.logn(data = data_RAMS.Pt1.Ex2,pp = "Blom",xlabel1 = "Cycles-to-Failure",stressunit1 = "MPa")
probplot.wbl(data = data_RAMS.Pt1.Ex2,pp = "Blom",xlabel1 = "Cycles-to-Failure",stressunit1 = "MPa")

# B. Maximum Likelihood EStimation
# Example 3: Run MLE using data from Example 1
probplot.logn(data_RAMS.Pt1.Ex1,"Blom","Time to Failure (hours)",MLE_i = 1, stressunit1 = "K")
# Then compare with the Weibull
probplot.wbl(data_RAMS.Pt1.Ex1,"Blom","Time to Failure (hours)",MLE_i = 1, stressunit1 = "K")

# C. Bayesian Updating
# Example 4.a: Update failure rate from uninformative prior
priorset <- c("uniform(0,1)")
Dat.out<-distribution.BAYESest(pt_est = c(0.05),
                               "Exponential",TTF= c(130,1400),
                               confid = 0.9,
                               priors = priorset,
                               nsamples = 30000,
                               burnin = 1000,nchains = 4)
# Example 4.b: Update failure rate from new prior (last posterior)
Post.dat <- extract(Dat.out$posterior.fit,c("lambda"))$lambda
distribution.fit(Post.dat)
priorset.2 <- c("gamma(2.964625, 0.0006598385)")
Dat.out.2<-distribution.BAYESest(pt_est = c(0.05),
                                 "Exponential",TTF= c(430,534,560,560,1403,2020),Tc = rep(2160,4),
                                 confid = 0.9,
                                 priors = priorset.2,
                                 nsamples = 20000,
                                 burnin = 1000,nchains = 4)

# ==============================================================================
# PART 2: Accelerate Life Testing
# ==============================================================================
# Example 1: ALT LSQ
lifestress.LSQest(data_RAMS.Pt1.Ex1,ls="Exponential2",dist = "Weibull",pp = "Blom",
                  xlabel1 = "Time to Failure (hours)",Suse = 353,
                  Llab = "Component Life (hours)",Slab = "Inverse Temperature (1/K)",stressunit1 = "K")

# Example 2: ALT MLE
lifestress.MLEest(data_RAMS.Pt1.Ex1,ls="Exponential2",dist = "Weibull",
                  Suse=353,xlabel1 = "Time to Failure (hours)",
                  Llab = "Component Life (hours)",Slab = "Inverse Temperature (1/K)",stressunit1 = "K")

# Example 3: Physics-of-Failure application - S-N Diagram and ALT MLE with Fatigue Test in Material Lab
data_RAMS.Pt2.Ex2 <- cbind(c(45000, 240000, 800000, 1500000, 2700000, 7800000, 10000000,26000000, 12000000, 22000000),
                           c(rep(1,7),rep(0,3)),
                           c(544.00, 510.35, 469.95, 436.23, 427.82, 411.00, 411.00,404.24, 397.48, 395.83))
probplot.logn(data_RAMS.Pt2.Ex2,"Blom","Cycles to Failure",stressunit1 = "MPa",MLE_i = 1)

lifestress.LSQest(data_RAMS.Pt2.Ex2,ls="InversePower2",dist = "Lognormal",pp = "Blom",
                  xlabel1 = "Cycles to Failure",Suse = 379,
                  Llab = "Cycles to Failure",Slab = "Alternating Stress (MPa)",stressunit1 = "MPa")
lifestress.MLEest(data_RAMS.Pt2.Ex2,ls="InversePower2",dist = "Lognormal",
                  Suse=379,xlabel1 = "Cycles to Failure",
                  Llab = "Cycles to Failure",Slab = "Alternating Stress (MPa)",stressunit1 = "MPa")

# Example 4: Physics-of-Failure application - S-N Diagram and ALT Bayesian with Accelerated Test data from components made from material
priorset <- c("lognormal(-0.715395,0.358001)","lognormal(2.9654,0.10560)","lognormal(-134.0451,8.199294)")
Bayes.Out<-lifestress.BAYESest(c(1,20,0.0001),ls="InversePower2",dist = "Lognormal",
                               TTF = c(2900000, 1400000, 9000000),SF = c(410.93, 404.72, 387.49),
                               Tc = c(10000000, 10000000, 10000000), Sc = c(394.38, 381.28, 310.95),
                               SUSE = 379, SACC = 400,confid = 0.9,
                               priors = priorset,
                               nsamples = 30000, burnin = 1000, nchains = 4)

# ==============================================================================
# PART 3: Accelerated Degradation Testing
# ==============================================================================
# ADT Example 1
data_RAMS.Pt3.Ex1 <- cbind(data.frame(rep(c(10, 20, 30, 40, 50, 60, 70),5),
                                      c(240.7, 316.2, 380.7, 409.1, 419.2, 455.8, 475.8,
                                        255.5, 297.9, 392.4, 436, 459.8, 472.8, 479.8,
                                        293.3, 316.1, 343.7, 397.5, 429.2, 444, 449.3,
                                        268.7, 305, 338.2, 359.5, 442.5, 464.6, 468.9,
                                        295.5, 348.1, 362, 400.9, 445.3, 455.3, 470.2),
                                      c(rep("Switch A",7),rep("Switch B",7),rep("Switch C",7),
                                        rep("Switch D",7),rep("Switch E",7)),
                                      rep(5,35)))
adt.rank(data_RAMS.Pt3.Ex1)

deg.output1 <- degradationlife.LSQest(data=data_RAMS.Pt3.Ex1,dl="SquareRoot2",dist="Normal",
                                      pp="Blom",D0=550,
                                      xlabel = "Time (Hours)",ylabel = "Resistance (Ohms, Ω)")
# May use either probability plotting (with MLE on) or by using the distribution.fit tool
distribution.fit(deg.output1$pseudo.time.to.failure.by.stress[,1])
probplot.exp2P(deg.output1[[2]],"Blom","Time to Failure (hours)",stressunit1 = "V",MLE_i = 1)

# ADT Example 2
data_RAMS.Pt3.Ex2 <- read.csv("https://raw.githubusercontent.com/Center-for-Risk-and-Reliability/RMT/main/CSVExampleData/Degradation_Data_1_Mass_Loss_by_Weight_gms_Example_5_2.csv")

data_RAMS.Pt3.Ex2 <- cbind(data.frame(rep(c(2,5,10,20,50,100,200,500),12),
                                      c(3.2, 4.1, 4.5, 4.7, 5.8, 6.8, 7.7, 9.6,
                                        2.7, 3.4, 3.8, 3.9, 5.4, 5.7, 6.3, 8.4,
                                        2.1, 2.7, 3.1, 3.3, 4, 4.6, 5.7, 6.6,
                                        2.6, 3.5, 4, 4, 5.2, 6.1, 6.7, 8.5,
                                        7.5, 7.8, 8.2, 10.6, 12.6, 13.3, 12.9, 14.8,
                                        7.5, 8.1, 9.8, 10.9, 14.8, 16.1, 17.3, 20.2,
                                        7, 8.9, 9.4, 11.1, 12.4, 13.5, 16.7, 17.3,
                                        7.8, 8.9, 10, 11.5, 13.7, 16.2, 16.2, 21,
                                        12.5, 15.4, 17.2, 20.5, 24.1, 27, 29.4, 37.9,
                                        11, 13.9, 16.1, 18.6, 22.2, 27.8, 31, 36.6,
                                        13, 15.1, 18.6, 20.2, 23.9, 29.7, 31.5, 39.6,
                                        11.7, 13.7, 16.7, 17.5, 22.3, 25.3, 32, 38.2),
                                      c(rep('Unit 01',8),rep('Unit 02',8),rep('Unit 03',8),rep('Unit 04',8),
                                        rep('Unit 05',8),rep('Unit 06',8),rep('Unit 07',8),rep('Unit 08',8),
                                        rep('Unit 09',8),rep('Unit 10',8),rep('Unit 11',8),rep('Unit 12',8)),
                                      c(rep(10,32),rep(50,32),rep(100,32))))
degradationlife.LSQest(data=data_RAMS.Pt3.Ex2,dl="Power",dist="Normal",
                       pp="Blom",D0=50, modelstress = "InversePower2",
                       xlabel = "Cycles to Failure (x100)",ylabel = "Wear (micro grams)")

# ADT Example 3 (MLE)
degradationlife.MLEest(data=data_RAMS.Pt3.Ex2,dl="Power",dist="Normal",pp="Blom",D0=50,
                       modelstress="InversePower2", confid=0.90,
                       xlabel="Cycles to Failure (x100)",ylabel="Wear (micro grams)",
                       stressunit1 = "grams")

degradationlife.MLEest(data=data_RAMS.Pt3.Ex2,dl="Power",dist="Normal",pp="Blom",D0=50,
                       modelstress="InversePower2", confid=0.90,Suse = 5,
                       xlabel="Cycles to Failure (x100)",ylabel="Wear (micro grams)",
                       stressunit1 = "grams")

# ADT Example 4 (Bayesian)
# Do this as a group because of its relative newness
# NEW DATA
data_RAMS.Pt3.Ex3 <- cbind(data.frame(c(rep(c(10,20,30,40,50,100,1000),4)),
                                      c(5.6, 5.7, 6.3, 6.7, 7, 10.2, 18.4,
                                        5.2, 5.5, 5.7, 6, 6.4, 9.3, 17.3,
                                        15, 15.6, 16.2, 16.8, 17.5, 22.3, 57.3,
                                        14.5, 14.7, 15.6, 17.4, 19.5, 21.5, 59.1),
                                      c(rep("Device A",7),rep("Device B",7),
                                        rep("Device C",7),rep("Device D",7)),
                                      c(rep(30,14),rep(70,14))))

priorset <- c("lognormal(-0.456903,0.1232223)",
              "normal(-0.661706,0.0140145)","lognormal(0.6961635,0.06053974)",
              "lognormal(-1.672606,0.01887561)","lognormal(-3.854416,0.08361823)")

Bayes.Out.2 <- degradationlife.BAYESest(pt_est = c(1,-1,1,1,1), data = data_RAMS.Pt3.Ex3,
                                        dl = "Power",dist="Normal", D0=50,
                                        modelstress="InversePower2",confid=0.95,
                                        SUSE=5,priors = priorset,
                                        nsamples = 30000, burnin = 1000, nchains = 4)


# ==============================================================================
# ==============================================================================
