# CHAPTER 4 PROBLEM 23
# Reuel Smith
# =================================================
data_SSALT_PROBLEM_23 <- cbind(c(2.01, 3.60, 4.1, 4.34,
                                 5.04, 5.94, 6.68, 7.09, 7.09, 7.49, 7.60, 8.23, 8.23, 8.23, 8.23, 12.05)*1000,
                               rep(1,16),
                               c(rep(12.18,4),rep(4.48,12)))
data_SSALT_PROBLEM_23a <- cbind(c(2.01, 3.60, 4.1, 4.34,rep(5,12),
                                 5.04, 5.94, 6.68, 7.09, 7.09, 7.49, 7.60, 8.23, 8.23, 8.23, 8.23, 12.05)*1000,
                               c(rep(1,4),rep(0,12),rep(1,12)),
                               c(rep(12.18,16),rep(4.48,12)))

table_SSALT_PROBLEM_23 <- cbind(c(12.18,4.48),c(5000,12500))
x0.P23<-stepstress.LSQest(data = data_SSALT_PROBLEM_23, stepstresstable = table_SSALT_PROBLEM_23,
                          ls = "InversePower2", dist = "Weibull", pp = "Blom", xlabel1 = "Cycles-to-Failure", stressunit1 = "MPa")
stepstress.MLEest(x0.P23$params_optimized, data = data_SSALT_PROBLEM_23, stepstresstable = table_SSALT_PROBLEM_23,
                  ls = "InversePower2", dist = "Weibull", confid = 0.9)
probplot.wbl(data_SSALT_PROBLEM_23a,pp = "Blom", xlabel1 = "Cycles-to-Failure (hours)",stressunit1 = "MPa")
