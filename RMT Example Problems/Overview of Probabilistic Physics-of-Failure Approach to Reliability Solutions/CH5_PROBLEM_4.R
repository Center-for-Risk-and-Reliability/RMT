# CHAPTER 5 PROBLEM 4 SOLUTIONS
time.P4 <- c(rep(1,13),rep(2,11),rep(3,9),rep(4,11))
deg.P4 <- c(437,446,497,503,705,737,748,788,818,860,875,934,1124,
            412,420,451,454,554,580,608,610,727,825,925,
            246,324,330,426,499,546,554,559,625,
            125,208,229,242,273,297,311,318,393,403,470)
params.P4  <- lm(deg.P4 ~ poly(time.P4, 1, raw=TRUE))
dlparams.P4 <- c(summary(params.P4)$coefficients[1,1],summary(params.P4)$coefficients[2,1])
R2.P4 <- summary(params.P4)$r.squared
SSE.P4 <- sum((fitted(params.P4) - deg.P4)^2)
sigma.P4 <- sqrt((1/(length(time.P4)-1))*SSE.P4)

# REGRESSION SCRIPT
CDF.1 <- 0.0019
CDF.2 <- 0.0028
CDF.3 <- 0.0065
CDF.4 <- 0.0583

wblparm.P4  <- lm(log(-log(1-c(CDF.1,CDF.2,CDF.3,CDF.4))) ~ poly(log(c(1,2,3,4)), 1, raw=TRUE))
beta.P4 <- summary(wblparm.P4)$coefficients[2,1]
alpha.P4 <- exp(-summary(wblparm.P4)$coefficients[1,1]/beta.P4)

