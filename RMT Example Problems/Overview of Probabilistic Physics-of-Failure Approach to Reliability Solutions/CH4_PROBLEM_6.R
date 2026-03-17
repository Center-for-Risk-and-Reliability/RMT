# CHAPTER 4 PROBLEM 6
# Gabriel Cruz and Reuel Smith
# =================================================
data_ALT_PROBLEM_6 <- cbind(c(248,456,500,500,500,500,
                              164,176,200,200,200,200,
                              164,176,180,180,180,180),
                            c(1,1,0,0,0,0,
                              1,1,0,0,0,0,
                              1,1,0,0,0,0),
                            c(406,406,406,406,406,406,
                              416,416,416,416,416,416,
                              426,426,426,426,426,426))
lifestress.MLEest(data = data_ALT_PROBLEM_6, ls = "Exponential2",dist = "Weibull",
                  confid = 0.95,pp = "Blom",xlabel1 = "Time (hours)",
                  Llab = "Life",Slab = "Temp",stressunit1 = "K")
