# Change Log for RMT
**Version 1.3.0.2** (September 19th, 2024)
* Bug fix update to the Stress-Strain Parameters tool (**stress_strain.params**). An issue causing several unnecessary warnings has been resolved.  Also the tool no longer calls nls.multstart library.
  
**Version 1.3.0.1** (September 18th, 2024)
* Bug fix update to the S-N diagram tool (**SN.diagram**). An issue with the **Ntrace** and **Strace** upper and lower estimates for equivalent alternating stress and equivalent fatigue cycle output respectively has been resolved.  Also an issue causing several unnecessary warnings has been resolved. 

**Version 1.3.0.0** (September 11th, 2024)
* Minor update to the S-N diagram tool (**SN.diagram**).  Tool now includes a confidence bound input (confid) which determines the upper and lower bounds of the S-N curve per the given data.  Tool also can calculate the upper and lower estimates of a given stress and/or fatigue cycle-to-failure trace.  That is, a fatigue cycle trace (**Ntrace**) will compute the equivalent alternating stress and its upper and lower confidence bounds and an alternating stress trace (**Strace**) will do likewise with the equivalent fatigue cycle and its upper and lower confidence bounds.
* Major update to probability plotting tools.  RMT can now estimate Three Parameter or Generalized Gamma Life Distribution (probplot.gam3P and probplotparam.gam3P) making eleven probability plot types total.  Some estimate refinements were performed on the Gamma and Three Parameter Weibull plotting and estimation tools.  Notably for the Gamma probability plot it is now a single probability plot.  Finally, all probability plots can be toggled between least square estimates and maximum likelihood estimates.
* Major update to the ALT least squared, maximum likelihood, and Bayesian tools (lifestress.LSQest, lifestress.MLEest, and lifestress.BAYESest).  Use stress has now been added as an input (Suse) in order to produce an estimate for the use or nominal life.  For the Bayesian ALT tool this allows for posterior use life estimation calculation and plotting.  In addition, an accelerated stress (SACC) is a new input for lifestress.BAYESest which allows calculation for posterior acceleration factor at SACC (documentation update pending).
* Major update to the ADT least squared and maximum likelihood tools (degradationlife.LSQest and degradationlife.MLEest).
* Minor updates to the step-stress ALT tool for least squares and maximum likelihood.  First, adt.full.LSQ and adt.full.MLE have been renamed degradationlife.LSQest and degradationlife.MLEest respectively (documentation update pending).  The degradationlife.LSQest tool plots the damage accumulation paths by default with some additional plot fixes planned for later this semester.  Both tools have a stress-model input (modelstress) which may be used to relate one of the damage-life model parameters to one of the RMT's available life-stress model.  Documentation to be updated by the next update.
* The Maximum Likelihood Variance-Covariance Matrix Select (MLE.var.covar.select) tool has a significant update to its MLE optimization tool selection process.  The tool now has a third input to indicate whether or not the log-likelihood has a bounded parameter or not.  This is necessary when dealing with the three-parameter Weibull and Generalized Gamma distributions because of certain boundary conditions (e.g. the location parameter of the three-parameter Weibull is a specific upper bound determined by the data).  For this, the three existing optimization tools (nlm, ucminf, and optim)  are set aside for general log-likelihoods.  Three new optimization tools have been added for bounded log-likelihoods: (1) nlminb, (2) DEoptim, and (3) nmkb (documentation update pending).

**Version 1.2.6.17** (December 6th, 2023)
* Made some updates to the life-stress and step-stress tools and made some changes to the help files for the adt.full.BAYES and stepstress.LSQest.
* Version 1.2.6.9 (November 15th, 2023)
* Made a bug fix to the Eyring and Eyring2 life stress models for  lifestress.BAYESest.

**Version 1.2.6.5-8** (November 7th, 2023)
* Added some missing scripts and updated help files for lifestress.LSQest, lifestress.MLEest, lifestress.BAYESest, stepstress.LSQest, stepstress.MLEest, and stepstress.BAYESest

**Version 1.2.6.4** (October 27th, 2023)
* Fixed an error in probplot.logn that affected its output and added fatigue stress concentration factor as an optional input ("options=list(Kf = )") in notch.lifestrain.life.trace.  Please see stress.concentration.factor help file for geometry definitions (also updated).

**Version 1.2.6.3** (October 23rd, 2023)
* Added missing distribution.MLEest tool for general MLE distribution evaluation.

**Version 1.2.6.2** (October 18th, 2023)
* Added reliability distribution MLE fitting tools for Weibull, Normal, Lognormal, Exponential, Logistic, Log-Logistic, Gumbel, and Gamma.
Updated output for Bayesian estimation tools distribution.BAYESest and lifestress.BAYESest.

**Version 1.2.6.1** (September 27th, 2023)
* Made a significant update to the stress.concentration.factor tool to include several new geometries and loading options (tension as well as bending and torsion).  Also updated the help file which now includes figures of geometries for reference.

**Version 1.2.6** (August 29th thru September 8th, 2023)
* Updated all probability plotting tools to include a column for non-parametric probabilistic density.
* Updated use of Bayesian based evaluations
* Added some new tools for some wear, creep, and corrosion applications:
> wear.sliding - Computes sliding wear output for a given scenario
> 
> wear.coefficient - Selects sliding wear coefficient for two materials from a database of values
> 
> creep.analysis - Performs creep data analysis by reading thermal acceleration creep data and determine creep model parameters (for a given creep model) and> creep rupture time with a given reference used for the model parameters
> 
> corrosion.pitting - Calculates pitting corrosion model parameters, rates, and times based on pit corrosion data.
> 
* Added some new tools for basic distribution fitting
* Added new probability plot options for Gamma, Gumbel, Logistic, and Log-logistic distributions.

**Version 1.2.5** (November 30th, 2022)
* Minor updates to adt.full.BAYES, adt.full.LSQ, adt.full.MLE, adt.full.multME.MLE, and adt.rank tools.

**Version 1.2.4** (November 2nd, 2022)
* Updated several help files for adherence to the  format that is now used in R documentation and I have added necessary references where needed.
* Also updated lognprobplot and stepstressLSQest for a continuity change.
* Please also note that there is a very important notice regarding Stan and R that may need to be addressed if you wish to use any of the Bayesian tools at the moment.

**Version 1.2.3** (October 19th, 2022)
* Corrected seven probability plot parameter generators to correct a bug caused by entering a single point of data or a series of data at the same value x.

**Version 1.2.2** (October 11th, 2022)
* Updated the six primary probability plotting tools: probplot.exp, probplot.exp2P, probplot.nor, probplot.logn, probplot.wbl, and probplot.wbl3P.  The plots now use ggplot2 for visualization allowing for data labeling and clearer reading of probability plots.

**Version 1.2.1** (October 10th, 2022)
* Corrected an output error on SN.diagram and updated the help files for 22 probability plotting-based tools.  Help files now updated to apply LaTeX scripting for equations.

**Version 1.2.0** (October 7th, 2022)
* Added seven new scripts:
> notch.SN.diagram - Generates the stress-life (or S-N) diagram based on various stress and fatigue life input data as well as notch effect.
>
> SN.meanstresseffect - Computes the equivalent stress amplitude  for a non-zero mean stress amplitude  or non-fully reversible scenario.
>
> stress.concentration.factor - Computes the stress concentration factor  for a specified geometry and set of dimensions.
>
> notch.lifestrain.life.trace - Computes the strain-life (in fatigue cycles  and reversals) for a given set of strain-life parameters and conditions with consideration for notch effect.
>
> fatigue.stress.concentration.factor.ratio - Approximates the fatigue notch factor ratio  based on the ultimate tensile strength (based on Figure 4.10 Bannantine).
>
> crack.propagation - Basic crack propagation analysis calculation.
>
> crack.correction - Computes the crack correction factor  or  for a specified geometry and set of dimensions.

**Version 1.1.0** (September 19th, 2022)
* Updated SN.diagram tool to run loading cases where loading is not fully reversed (that is ) and combined loading conditions on 3D or multiaxial fatigue cases.
* Added four new scripts related to stress-life and strain-life evaluation (stress_strain.params, hysteresisloop.plot, lifestrain.life.trace) and Miner's Rule (var.amp.loadingdamage.model).

**Version 1.0.0** (September 1st, 2022)
* Initial Reliability Modeling Toolkit
