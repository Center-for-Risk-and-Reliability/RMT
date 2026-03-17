# CHAPTER 4 PROBLEM 14
# Gabriel Cruz and Reuel Smith
# =================================================
data_ALT_PROBLEM_14 <- cbind(c(780,805,900,1020,
                 490,530,580,615,
                 290,310,355,390,
                 180,200,250),
               c(1,1,1,1,
                 1,1,1,1,
                 1,1,1,1,
                 1,1,1),
               c(340,340,340,340,
                 340,340,340,340,
                 370,370,370,370,
                 370,370,370),
               c(4,4,4,4,
                 6,6,6,6,
                 4,4,4,4,
                 6,6,6))
lifestress.LSQest(data = data_ALT_PROBLEM_14, ls = "TempNonthermal", dist = "Lognormal", pp = "Blom",
                  stressunit1 = "K", stressunit2 = "kN")
lifestress.MLEest(data = data_ALT_PROBLEM_14, ls = "TempNonthermal", dist = "Lognormal", pp = "Blom", confid = 0.9,
                  stressunit1 = "K", stressunit2 = "kN")
