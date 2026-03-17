# CHAPTER 4 PROBLEM 18
# Gabriel Cruz
# =================================================
# 50 units tested at stress levels 150, 175, 200, 250, and 300 °C

data_ALT_PROBLEM_18 <- cbind(c(788,rep(1536,3),rep(2304,5),rep(2304,41),
                               rep(384,4),rep(788,27),rep(1536,16),rep(1536,3)),
                             c(rep(1,9),rep(0,41),
                               rep(1,47),rep(0,3)),
                             c(rep(250+273.15,50),
                               rep(300+273.15,50)))
# Part a)
probplot.wbl(data = data_ALT_PROBLEM_18,pp = "Blom",
             xlabel1 = "Time to Failure (hours)",stressunit1 = "K")
lifestress.LSQest(data = data_ALT_PROBLEM_18, ls = "Arrhenius", dist = "Weibull",pp = "Blom",
                  xlabel1 = "Time to Failure (hours)",Llab = "Life (hours)", Slab = "Inverse Temperature (1/K)",
                  stressunit1 = "K", Suse = 273.15+40)

# Part b)
probplot.wbl(data = data_ALT_PROBLEM_18,pp = "Blom",
             xlabel1 = "Time to Failure (hours)",MLE_i = 1,stressunit1 = "K")
lifestress.MLEest(data = data_ALT_PROBLEM_18, ls = "Arrhenius", dist = "Weibull",pp = "Blom", confid = 0.9,
                  xlabel1 = "Time to Failure (hours)",Llab = "Life (hours)", Slab = "Inverse Temperature (1/K)",
                  stressunit1 = "K", Suse = 273.15+40)
# Part e)
out.P18 <- lifestress.BAYESest(pt_est = c(2,0.1,0.01), ls = "Arrhenius", dist = "Weibull",
                    TTF = c(788,rep(1536,3),rep(2304,5),rep(384,4),rep(788,27),rep(1536,16)),
                    SF = c(rep(250+273.15,9), rep(300+273.15,47)),
                    Tc = c(rep(2304,41),rep(1536,3)),
                    Sc = c(rep(250+273.15,41), rep(300+273.15,3)), SUSE = 313.15,
                    priors = c("uniform(1,5)","uniform(0.001,1)","uniform(0.00001,1)"),
                    nsamples = 20000,burnin = 1000,nchains = 4)
