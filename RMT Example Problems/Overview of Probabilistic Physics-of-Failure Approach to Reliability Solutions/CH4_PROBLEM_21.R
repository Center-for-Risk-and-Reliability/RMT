# CHAPTER 4 PROBLEM 21
# Reuel Smith
# =================================================
data_SSALT_PROBLEM_21 <- cbind(c(252,280,
                                 320,328,335,
                                 354,361,362,368,375,375,375),
                               c(rep(1,2),
                                 1,0,1,
                                 1,1,1,1,rep(0,3)),
                               c(rep(175,2),
                                 rep(200,3),
                                 rep(250,7)))
table_SSALT_PROBLEM_21 <- cbind(c(125,175,200,250),
                                c(200,100,50,25))
x0.P21<-stepstress.LSQest(data = data_SSALT_PROBLEM_21, stepstresstable = table_SSALT_PROBLEM_21,
                          ls = "InversePower2", dist = "Weibull", pp = "Blom", xlabel1 = "Time-to-Failure (hours)",
                          stressunit1 = "%")
probplot.wbl(x0.P21[[5]],pp = "Blom", xlabel1 = "Time-to-Failure (hours)",stressunit1 = "%")

