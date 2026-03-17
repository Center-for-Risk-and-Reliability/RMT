# CHAPTER 4 PROBLEM 9
# Reuel Smith
# =================================================
data_SSALT_PROBLEM_9 <- cbind(c(0.56e6,1.93e6,
                              2.12e6,2.81e6,
                              3.21e6,3.37e6,3.72e6,3.91e6,
                              rep(4e6,7)),
                            c(rep(1,8),rep(0,7)),
                            c(rep(120,2),rep(120*1.1,2),rep(120*1.2,11)))
table_SSALT_PROBLEM_9 <- cbind(c(120,120*1.1,120*1.2),c(2e6,1e6,1e6))
stepstress.LSQest(data = data_SSALT_PROBLEM_9, stepstresstable = table_SSALT_PROBLEM_9,
                  ls = "InversePower2", dist = "Weibull", pp = "Blom", xlabel1 = "Cycles-to-Failure", stressunit1 = "MPa")

dat_plot_0 <- stepstress.LSQest(data = data_SSALT_PROBLEM_9, stepstresstable = table_SSALT_PROBLEM_9,
                  ls = "InversePower2", dist = "Weibull", pp = "Blom", xlabel1 = "Cycles-to-Failure", stressunit1 = "MPa")[[5]]
probplot.wbl(dat_plot_0,pp = "Blom", xlabel1 = "Cycles-to-Failure",stressunit1 = "MPa")

lifestress.LSQest(data = dat_plot_0, ls = "InversePower2", dist = "Weibull", pp = "Blom",
                  xlabel1 = "Cycles-to-Failure", Llab = "Characteristic Life (cycles)", Slab = "Alternating Stress S_a (MPa)",
                  stressunit1 = "MPa")

log(8.653722e-39)
1/14.5995
(exp(-87.64283))^-0.0684955

(120/132)^(-1/0.0684955)
(120/144)^(-1/0.0684955)
(132/144)^(-1/0.0684955)

2e6/4.0208
2e6/14.32209
1e6/3.562001

1e6 + 497413.5
1e6 + 280741.1 + 139644.4

AdjAF.1.to.2 <- (27/13)*(2e6/(1e6 + 497413.5))
AdjAF.1.to.3 <- (53/13)*(2e6/(1e6 + 280741.1 + 139644.4))
AdjAF.2.to.3 <- (53/27)*((1e6 + 497413.5)/(1e6 + 280741.1 + 139644.4))

1/(log(2.774014)/log(120/132))
1/(log(5.740587)/log(120/144))
1/(log(2.069415)/log(132/144))

A <- exp(log(133) + 0.1057951*log(1420386) - (0.1057951/0.8522432)*log(-log(.5)))
