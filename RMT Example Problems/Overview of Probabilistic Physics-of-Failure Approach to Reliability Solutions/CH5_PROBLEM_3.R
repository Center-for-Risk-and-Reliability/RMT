# CHAPTER 5 PROBLEM 3
# Reuel Smith
# =================================================
data_ADT_PROBLEM_3 <- cbind(data.frame(time=rep(c(0,50,100,150,200,250),12),
                            color_fading=c(0.16, 0.21, 0.27, 0.38, 0.55, 0.71,
                                        0.05, 0.12, 0.18, 0.3, 0.46, 0.62,
                                        0.11, 0.15, 0.21, 0.33, 0.49, 0.66,
                                        0.12, 0.17, 0.23, 0.34, 0.5, 0.67,
                                        0.08, 0.33, 0.45, 0.8, 1.4, 2.15,
                                        0.12, 0.25, 0.54, 0.95, 1.59, 2.22,
                                        0.13, 0.18, 0.53, 1.06, 1.89, 2.4,
                                        0.08, 0.27, 0.64, 1.23, 1.97, 2.36,
                                        0.06, 0.59, 1.15, 1.85, 2.3, 2.7,
                                        0.12, 0.65, 1.35, 1.92, 2.35, 2.85,
                                        0.14, 0.62, 1.29, 1.87, 2.22, 2.74,
                                        0.09, 0.52, 1.32, 1.97, 2.45, 2.97),
                            unit_no=c(rep('Unit 1',6),rep('Unit 2',6),rep('Unit 3',6),rep('Unit 4',6),
                                      rep('Unit 5',6),rep('Unit 6',6),rep('Unit 7',6),rep('Unit 8',6),
                                      rep('Unit 9',6),rep('Unit 10',6),rep('Unit 11',6),rep('Unit 12',6)),
                            UV_Irradiance=c(rep(1,24),rep(3,24),rep(5,24))))
adt.rank(data_ADT_PROBLEM_3)  # Check ranking of data against other models
# $adt.rank.average
#              Average Rank
# Linear                  1
# Exponential             4
# Square-Root             2
# Square-Root2            3
# Power                   5
# Logarithmic             6
# Lloyd-Lipow             7
# Mitsuom                 8

degradationlife.LSQest(data=data_ADT_PROBLEM_3,dl="Linear",dist="Normal",
                       pp="Blom",Df=1, D0 = NULL, modelstress = NULL,
                       xlabel = "Time (hours)",ylabel = "Color Fading (ΔE)",stressunit1 = "W/m²")
degradationlife.LSQest(data=data_ADT_PROBLEM_3,dl="SquareRoot",dist="Normal",
                       pp="Blom",Df=1, D0 = NULL, modelstress = NULL,
                       xlabel = "Time (hours)",ylabel = "Color Fading (ΔE)",stressunit1 = "W/m²")
# ===================================================================
# We go with the Square Root degradation life model  D(l) = (a + bl)²
# ===================================================================
# Now check for best modelstress relation
degradationlife.LSQest(data=data_ADT_PROBLEM_3,dl="SquareRoot",dist="Normal",
                       pp="Blom",Df=1, D0 = NULL, modelstress = "Linear",
                       xlabel = "Time (hours)",ylabel = "Color Fading (ΔE)",stressunit1 = "W/m²")
degradationlife.LSQest(data=data_ADT_PROBLEM_3,dl="SquareRoot",dist="Normal",
                       pp="Blom",Df=1, D0 = NULL, modelstress = "InversePower2",
                       xlabel = "Time (hours)",ylabel = "Color Fading (ΔE)",stressunit1 = "W/m²")
degradationlife.LSQest(data=data_ADT_PROBLEM_3,dl="SquareRoot",dist="Normal",
                       pp="Blom",Df=1, D0 = NULL, modelstress = "Exponential",
                       xlabel = "Time (hours)",ylabel = "Color Fading (ΔE)",stressunit1 = "W/m²")
# ===================================================================
# Of linear a(S) = (b_0 + S*a_0), inverse power a(S) = 1/[b_0*(S^a_0)],
# and exponential a(S) = b_0*exp(a_0*S), we go with the exponential
# parameter-stress relation because it has the lowest SSE
# D(l,S) = (b_0*exp(a_0*S) + bl)²
# ===================================================================
degradationlife.MLEest(data=data_ADT_PROBLEM_3,dl="SquareRoot",dist="Normal",
                       pp="Blom",Df=1, D0 = NULL, modelstress = "Exponential",confid = 0.90,Suse = 0.25,
                       xlabel = "Time (hours)",ylabel = "Color Fading (ΔE)",stressunit1 = "W/m²")

