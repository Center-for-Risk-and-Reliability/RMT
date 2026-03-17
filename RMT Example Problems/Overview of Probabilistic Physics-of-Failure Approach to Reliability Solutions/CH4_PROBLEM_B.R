# CHAPTER 4 PROBLEM Possible B
# Gabriel Cruz and Reuel Smith
# =================================================
data_ALT_PROBLEM_11 <- cbind(c(1000,1200,1500,rep(1500,7),900,1300,1400,rep(1400,7),500,800,1100,rep(1100,7),
                               700,900,1400,rep(1400,7),500,700,1100,rep(1100,7),300,800,1000,rep(1000,7),
                               500,600,1100,rep(1100,7),300,400,1000,rep(1000,7),500,600,900,rep(900,7)),
                             c(1,1,1,rep(0,7),1,1,1,rep(0,7),1,1,1,rep(0,7),
                               1,1,1,rep(0,7),1,1,1,rep(0,7),1,1,1,rep(0,7),
                               1,1,1,rep(0,7),1,1,1,rep(0,7),1,1,1,rep(0,7)),
                             c(rep(303.15,10),rep(313.15,10),rep(323.15,10),
                               rep(303.15,10),rep(313.15,10),rep(323.15,10),
                               rep(303.15,10),rep(313.15,10),rep(323.15,10)),
                             c(rep(100,30),rep(200,30),rep(300,30)))
lifestress.LSQest(data = data_ALT_PROBLEM_11, ls = "TempNonthermal", dist = "Weibull", pp = "Blom",
                  stressunit1 = "K", stressunit2 = "V")
lifestress.MLEest(data = data_ALT_PROBLEM_11, ls = "TempNonthermal", dist = "Weibull", pp = "Blom", confid = 0.95,
                  stressunit1 = "K", stressunit2 = "V")
probplot.wbl(data = data_ALT_PROBLEM_11, pp = "Blom", xlabel1 = "Time to Failure (hours)", MLE_i = 1,
             stressunit1 = "K", stressunit2 = "V")
