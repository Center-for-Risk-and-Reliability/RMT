# CHAPTER 4 PROBLEM 22
# Reuel Smith
# =================================================
data_SSALT_PROBLEM_22 <- cbind(c(380,
                                 480,590,
                                 677, 710, 735,
                                 788, 810, 853, 888,
                                 912, 923, 940, 948, 955, 972),
                               rep(1,16),
                               c(0.40,rep(0.60,2),rep(0.70,3),rep(0.80,4),rep(0.90,6)),
                               c(0.1,rep(0.2,2),rep(0.3,3),rep(0.4,4),rep(0.45,6)))
data_SSALT_PROBLEM_22.fit <- cbind(c(380,
                                 80,190,
                                 77, 110, 135,
                                 38, 69, 103, 138,
                                 12, 23, 40, 48, 55, 72),
                               rep(1,16),
                               c(0.40,rep(0.60,2),rep(0.70,3),rep(0.80,4),rep(0.90,6)),
                               c(0.1,rep(0.2,2),rep(0.3,3),rep(0.4,4),rep(0.45,6)))
table_SSALT_PROBLEM_22 <- cbind(c(0.40,0.60,0.70,0.80,0.90),c(0.1,0.2,0.3,0.4,0.45),c(400,200,150,150,100))
probplot.wbl(data = data_SSALT_PROBLEM_22.fit,pp = "Blom",xlabel1 = "Time-to-Failure (hours)", stressunit1 = "RH", stressunit2 = "g", MLE_i = 1)
probplot.logn(data = data_SSALT_PROBLEM_22.fit,pp = "Blom",xlabel1 = "Time-to-Failure (hours)", stressunit1 = "RH", stressunit2 = "g", MLE_i = 1)
probplot.gam(data = data_SSALT_PROBLEM_22.fit,pp = "Blom",xlabel1 = "Time-to-Failure (hours)", stressunit1 = "RH", stressunit2 = "g", MLE_i = 1)
probplot.logist(data = data_SSALT_PROBLEM_22.fit,pp = "Blom",xlabel1 = "Time-to-Failure (hours)", stressunit1 = "RH", stressunit2 = "g", MLE_i = 1)

stepstress.LSQest(data = data_SSALT_PROBLEM_22, stepstresstable = table_SSALT_PROBLEM_22,
                  ls = "TempNonthermal", dist = "Weibull", pp = "Blom", xlabel1 = "Time-to-Failure (hours)", stressunit1 = "RH", stressunit2 = "g")

lifestress.LSQest(data = gogo, ls = "TempNontherm", dist = "Weibull",pp = "Blom",xlabel1 = "Time to Failure (hours",stressunit1 = "RH", stressunit2 = "g")

39.1345932*(0.1^-1.5324147)*exp(0.4318075/0.4)*((-log(1-0.5))^(1/1.4840672))
