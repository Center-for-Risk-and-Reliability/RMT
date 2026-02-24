# CHAPTER 5 PROBLEM 2
# Reuel Smith
# =================================================
data_ADT_PROBLEM_2 <- cbind(data.frame(time=c(rep(c(0,1,2,3,4,5,6),3)),
                                      fading=c(0, -0.05, -0.15, -0.22, -0.31, -0.35, -0.41,
                                               0, -0.03, -0.11, -0.18, -0.28, -0.34, -0.43,
                                               0, -0.01, -0.09, -0.19, -0.29, -0.35, -0.39),
                                      sample=c(rep("Dye 1 - Magenta",7),rep("Dye 2 - Yellow",7),
                                               rep("Dye 3 - Cyan",7)),
                                      temp=c(rep(297,21))))

degradationlife.LSQest(data = data_ADT_PROBLEM_2, dl = "Linear", dist = "Normal",
                       pp = "Blom", Df = -0.5, D0 = 0, modelstress = NULL,
                       xlabel = "Time to Failure (years)",ylabel = "Fading (percent)")

degradationlife.MLEest(data = data_ADT_PROBLEM_2, dl = "Linear", dist = "Normal",
                       pp = "Blom",Df = -0.5, D0 = 0, modelstress = NULL, confid=0.90, Suse = NULL,
                       xlabel = "Time to Failure (years)",ylabel = "Fading (percent)")
