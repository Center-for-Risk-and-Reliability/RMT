# CHAPTER 4 PROBLEM 7
# Gabriel Cruz and Reuel Smith
# =================================================
data_ALT_PROBLEM_7 <- cbind(c(78,88,99,105,117,138,200,
                              660,711,798,911,1000,1012,
                              152,193,203,223,272,
                              311,380,417,518,611),
               c(1,1,1,1,1,1,0,
                 1,1,1,1,0,1,
                 1,1,1,1,1,
                 1,1,1,1,1),
               c(90,90,90,90,90,90,90,
                 30,30,30,30,30,30,
                 60,60,60,60,60,
                 50,50,50,50,50))
probplot.wbl(data = data_ALT_PROBLEM_7, pp = "Blom", xlabel1 = "Time to Failure (hours)", stressunit1 = "N")
lifestress.LSQest(data = data_ALT_PROBLEM_7, ls = "Linear", dist = "Weibull", pp = "Blom",
                  xlabel1 = "Time to Failure (hours)", Llab = "Characteristic Life (hours)", Slab = "Load, W (Newtons)",
                  stressunit1 = "N", Suse = 12)
