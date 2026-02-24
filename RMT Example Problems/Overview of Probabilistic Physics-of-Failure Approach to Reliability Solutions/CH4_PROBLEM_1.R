# CHAPTER 5 PROBLEM 2
# Gabriel and Reuel Smith
# =================================================
# Part a)
data_ALT_PROBLEM_1 <- cbind(c(300,340,345,349,361,362,363,369,374,379,380,390,
                             250,251,252,260,262,263,264,270,271,272,280,285,
                             200,205,207,209,210,211,215,220,222,225,228,230),
                           rep(1,36),
                           c(rep(200+273.15,12),
                             rep(250+273.15,12),
                             rep(200+273.15,12)),
                           c(rep(10,12),
                             rep(10,12),
                             rep(15,12)))    # Enter the data
# Run probability plotting tool probplot.wbl to get least squares estimates for each stress level
probplot.wbl(data = data_ALT_PROBLEM_1,pp = "Blom",
             xlabel1 = "Failure Time (Hours)",
             stressunit1 = "K", stressunit2 = "Intensity")

# Part 2)
# Now run the lognormal probability plotting tool probplot.logn
probplot.logn(data = data_ALT_PROBLEM_1,pp = "Blom",
             xlabel1 = "Failure Time (Hours)",
             stressunit1 = "K", stressunit2 = "Intensity")
