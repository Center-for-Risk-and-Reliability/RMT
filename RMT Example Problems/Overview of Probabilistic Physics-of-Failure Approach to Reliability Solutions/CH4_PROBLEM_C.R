# CHAPTER 4 PROBLEM Possible C
# Gabriel Cruz and Reuel Smith
# =================================================
data_ALT_PROBLEM_12 <- cbind(c(29.1,80.7,100,100,100,100,100,100,100,100,
                               47.5,71.8,73.7,86.2,100,100,100,100,100,100,
                               29.5,52,56.3,63.5,92.5,99.5,100,100,100,100,
                               26.1,47.5,53.4,56.1,61.8,66.6,76.6,80.9,100,100),
                             c(1,1,0,0,0,0,0,0,0,0,
                               1,1,1,1,0,0,0,0,0,0,
                               1,1,1,1,1,1,0,0,0,0,
                               1,1,1,1,1,1,1,1,0,0),
                             c(300,300,300,300,300,300,300,300,300,300,
                               350,350,350,350,350,350,350,350,350,350,
                               400,400,400,400,400,400,400,400,400,400,
                               500,500,500,500,500,500,500,500,500,500))

lifestress.MLEest(data = data_ALT_PROBLEM_12,ls = "Arrhenius",dist = "Lognormal",confid = 0.95,
                  xlabel1 = "Failure Time (hours)", stressunit1 = "K")
probplot.logn(data = data_ALT_PROBLEM_12, pp = "Blom", xlabel1 = "Failure Time (hours)", MLE_i = 1,
             stressunit1 = "K")
