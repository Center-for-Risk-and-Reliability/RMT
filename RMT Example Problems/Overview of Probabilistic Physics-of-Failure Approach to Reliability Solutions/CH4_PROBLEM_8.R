# CHAPTER 4 PROBLEM 8
# Gabriel Cruz and Reuel Smith
# =================================================
data_ALT_PROBLEM_8 <- cbind(c(3450,4340,4760,5320,5740,6160,7200,
                              3300,3720,4180,4560,4920,5280,6000,
                              2645,3100,3400,3800,4100,4400,4700),
                            c(1,1,1,1,1,1,0,
                              1,1,1,1,1,1,0,
                              1,1,1,1,1,1,0),
                            c(2.71,2.71,2.71,2.71,2.71,2.71,2.71,
                              2.81,2.81,2.81,2.81,2.81,2.81,2.81,
                              2.91,2.91,2.91,2.91,2.91,2.91,2.91))
lifestress.LSQest(data = data_ALT_PROBLEM_8, ls = "InversePower2", dist = "Weibull", pp = "Blom",
                  xlabel1 = "Time to Failure (hours)", Llab = "Characteristic Life (hours)", Slab = "Load, W (MPa)",
                  stressunit1 = "MPa")

lifestress.MLEest(data = data_ALT_PROBLEM_8, ls = "InversePower2", dist = "Weibull", pp = "Blom", confid = 0.95,
                  xlabel1 = "Time to Failure (hours)", Llab = "Characteristic Life (hours)", Slab = "Load, W (MPa)",
                  stressunit1 = "MPa")
