# CHAPTER 4 PROBLEM 5
# Gabriel Cruz and Reuel Smith
# =================================================
# Part b)
data_ALT_PROBLEM_5 <- cbind(c(2411, 2593, 3001,rep(3200,7),
                              800, 1231, 1231, 1317, 1362, 1362, 1371, 1431, 1771, 1966,
                              219, 248, 283, 335, 361, 380, 452, 452, 518, 611),
                            c(rep(1,3),rep(0,7),
                              rep(1,10),
                              rep(1,10)),
                            c(rep(150+273.15,10),
                              rep(175+273.15,10),
                              rep(200+273.15,10)))
probplot.logn(data = data_ALT_PROBLEM_5,pp = "Blom",xlabel1 = "Time (hours)",stressunit1 = "K")

# Part c)
probplot.wbl(data = data_ALT_PROBLEM_5,pp = "Blom",xlabel1 = "Time (hours)",stressunit1 = "K")

# Part f)
lifestress.LSQest(data = data_ALT_PROBLEM_5,ls = "InversePower2",dist = "Lognormal",pp = "Blom",
                  xlabel1 = "Failure Time (Hours)",
                  Llab = "Characteristic Life (Hours)", Slab = "Temperature (K)",
                  Suse = 80+273.15, stressunit1 = "K")

# Part g)
param.out <- lifestress.LSQest(data = data_ALT_PROBLEM_5,ls = "InversePower2",dist = "Lognormal",pp = "Blom",
                  xlabel1 = "Failure Time (Hours)",
                  Llab = "Characteristic Life (Hours)", Slab = "Temperature (K)",
                  Suse = 80+273.15, stressunit1 = "K")$LSQ.point.estimates
accelfactor("InversePower2",param.out[2:3],Suse = 353.15,Sacc = c(423.150,  448.150, 473.150))

