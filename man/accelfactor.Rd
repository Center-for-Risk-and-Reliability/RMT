\name{accelfactor}
\alias{accelfactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Acceleration Factor Calculator
}
\description{
Computes the acceleration factors for a given set of output parameters obtained by least-squares, maximum likelihood, or Bayesian posterior estimate.
}
\usage{
accelfactor(params, lsm, Sacc, Suse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{params}{
Full line of life-stress distribution parameters.  The output from any of the life-stress esimator functions (\code{\link{lifestress.LSQest}}, \code{\link{lifestress.MLEest}}, or \code{\link{lifestress.BAYESest}}).
}
  \item{lsm}{
The life-stress model definition for analysis.  Choose from the list of life-models outlined in \code{\link{lifestress.select}}.
}
  \item{Sacc}{
The accelerated stress (or stresses) for the accletated life test scenario.
  \itemize{
  \item If a single stress, \code{Sacc} is a single value input for one accelerated stress or a vector of accelerated stresses.
  \item If dual or multiple stress, enter as a vector starting with temperature based stress for one accelerated stress.  For multiple accelerated stress levels, enter each stress type as a vector of equal length and use cbind to create a matrix.
  }
}
  \item{Suse}{
The use stress for the accletated life test scenario.  If dual or multiple stress, enter as a vector starting with temperature based stress.
}
}
\details{
Acceleration factor \eqn{AF} calculations are based on the given stress-life equations \eqn{L(S)} and the standard notation where,

\deqn{AF = \frac{L(S_{use})}{L(S_{acc})}}
}
\value{
Returns the acceleration factor for this accletated life test scenario.
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2024)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: USe the failure data from Example 4.1 of the textbook to find tha acceleration factor at 400 Kelvin given the use temperature of 333 Kelvin.
time <- c(248,456,528,731,813,965,972,1528,164,176,289,319,386,459,92,105,155,184,219,235)
cens <- rep(1,20)
temp <- c(rep(406,8),rep(436,6),rep(466,6))
mat1 <- cbind(time,cens,temp)

# Estimate the LSQ parameters for MLE estimate
LSQparams<-lifestress.LSQest(mat1,"Arrhenius","Weibull","Blom")[[3]]
MLEparams<-lifestress.MLEest(LSQparams,"Arrhenius","Weibull",time,temp)[[1]]

# Note: Only need the last two estimates (the stress-life for the acceleration factor)
accelfactor(MLEparams[2:3],"Arrhenius",Sacc = 400, Suse = 333)

## Example 2: Find the acceleration factors (from the Example 1 scenario) at the test levels (406, 436, and 466 Kelvin)
accelfactor(MLEparams[2:3],"Arrhenius",Sacc = c(406,436,466), Suse = 333)

## Example 3: Use the 12 capacitor failure data from Example 4.11 of the textbook to find the acceleration factors at the three given stress levels (LEVEL 1 - 348 Kelvin/3 Volts, LEVEL 2 - 348 Kelvin/5 Volts, LEVEL 3 - 378 Kelvin/3 Volts) given the use stress of 333 Kelvin/2 Volts.  Assume a Lognormal life distribution and a temperature-nonthermal life-stress relation.
time2 <- c(620,632,658,822,380,416,460,596,216,246,332,400)
cens2 <- rep(1,12)
temp2 <- c(rep(348,8),rep(378,4))
volt2 <- c(rep(3,4),rep(5,4),rep(3,4))
mat2 <- cbind(time2,cens2,temp2,volt2)

LSQparams2<-lifestress.LSQest(mat2,"TempNonthermal","Lognormal","Blom")[[3]]
MLEparams2<-lifestress.MLEest(LSQparams2,"TempNonthermal","Lognormal",time2,cbind(temp2,volt2))[[1]]

# Note: Only need the last three estimates (the stress-life for the acceleration factor)
accelfactor(MLEparams2[2:4],"TempNonthermal",Sacc = cbind(c(348,348,378),c(3,5,3)), Suse = c(333,2))


## Example 4: Find the acceleration factor (from the Example 2 scenario) at 406 Kelvin/6 Volts.
accelfactor(MLEparams2[2:4],"TempNonthermal",Sacc = c(406,6), Suse = c(333,2))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
