# CHAPTER 4 PROBLEM 19
# Reuel Smith
# =================================================
data_SSALT_PROBLEM_19 <- cbind(c(289,310,320,
                                 362,368,369,
                                 371,372,378,380,380),
                               c(rep(1,9),rep(0,2)),
                               c(rep(5,3),rep(7,3),rep(9,5)))
table_SSALT_PROBLEM_19 <- cbind(c(2,5,7,9),c(250,100,20,10))
x0.P19<-stepstress.LSQest(data = data_SSALT_PROBLEM_19, stepstresstable = table_SSALT_PROBLEM_19,
                  ls = "InversePower2", dist = "Weibull", pp = "Blom", xlabel1 = "Time-to-Failure (hours)", stressunit1 = "V",Suse = 2)
probplot.wbl(x0.P19[[5]],pp = "Blom", xlabel1 = "Time-to-Failure (hours)",stressunit1 = "V")
# Part b
stepstress.MLEest(c(1.575039e+00, 4.240255e+00, 4.757519e-06),data = data_SSALT_PROBLEM_19, stepstresstable = table_SSALT_PROBLEM_19,
                  ls = "InversePower2", dist = "Weibull", Suse = 2, confid = 0.9)

# Part d
Bayes.out.P19 <- stepstress.BAYESest(c(1.5, 3, 0.000003162278),data = data_SSALT_PROBLEM_19, stepstresstable = table_SSALT_PROBLEM_19,
                    ls = "InversePower2", dist = "Weibull", confid = 0.95,SUSE = 2,
                    priors = c("uniform(0.5,3)","uniform(2,5)","uniform(0.0000003162278,0.00003162278)"),
                    nsamples = 20000,burnin = 1000,nchains = 4)
