\name{creep.analysis}
\alias{creep.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creep Analyzer
}
\description{
Performs creep data analysis by reading thermal acceleration creep data and determine creep model parameters (for a given creep model) and creep rupture time with a given reference used for the model parameters.
}
\usage{
creep.analysis(data, model, creepproperties, units)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
      Tabular form of creep test data where column 1 is temperature (in Celsius or Fahrenheit), column 2 is stress, and column 3 is the rupture time (in hours).
}
  \item{model}{
      Defines the creep model being fit to the data. Enter \code{1} for \emph{Larson-Miller}, \code{2} for \emph{Manson-Haferd}, and \code{3} for \emph{Sherby-Dorn}.
      \enumerate{
          \item \strong{Larson-Miller}:
          \deqn{P = \left(T+273.15\right)\left(\log_{10}t+C\right)\text{  (Metric or S.I.)}}
          \deqn{P = \left(T+459.67\right)\left(\log_{10}t+C\right)\text{  (English)}}

          where \eqn{P} is the Larson-Miller parameter, \eqn{T} is temperature in Celsius (Metric or S.I.) or Fahrenheit (English), and \eqn{C} is the material constant.

          \item \strong{Manson-Haferd}:
          \deqn{P = \frac{T-T_a}{\log_{10}t - \log_{10}t_a}}

          where \eqn{P} is the Manson-Haferd parameter, \eqn{T} is temperature in in Celsius (Metric or S.I.) or Fahrenheit (English), \eqn{T_a} is a material constant in Celsius (Metric or S.I.) or Fahrenheit (English), and \eqn{\log_{10}t} is the material constant.

          \item \strong{Sherby-Dorn}:
          \deqn{P = \log_{10}t - 0.43\frac{E_a}{k_B T}}

          or

          \deqn{P = log_{10}t - 0.43\frac{Q}{R T}}

          where \eqn{P} is the Sherby-Dorn parameter, \eqn{T} is temperature in Kelvin, \eqn{E_a} is the activation energy in eV (\eqn{Q} is activation energy in J/mol), \eqn{k_B} is the Boltzmann Constant \eqn{8.617\times 10^{-5}}, \eqn{R} is the universal Gas Constant \eqn{8.314} J/mol-K, and \eqn{C} is the material constant.
    }
}
  \item{creepproperties}{
      The material properties based on the creep model (see \strong{Details}), reference creep test temperature \code{Tref}, and/or stress trace \code{Strace} entered as a list.
}
  \item{units}{
      Defines the units for the measurements: \code{1} for Metric or SI (temperature \emph{'Celsius'}) and \code{2} for English (temperature \emph{'Fahrenheit'}).  Metric will be used as default if this heading is left empty.
}
}
\details{
Creep analysis requires adherence to one of several models of which this tool has three: \emph{Larson-Miller}, \emph{Manson-Haferd}, and \emph{Sherby-Dorn}.  Enter material properies in \code{creepproperties} list as follows:
  \itemize{
    \item \strong{Larson-Miller} - \code{creepproperties = list(C = ...)}
    \item \strong{Manson-Haferd} - \code{creepproperties = list(temp_a = ..., log10t_a = ...)}
    \item \strong{Sherby-Dorn} - \code{creepproperties = list(Q = ...)} or \code{list(Ea = ...)}
  }
}
\value{
Returns the following
 \item{CreepPlot}{
    A plot of the creep data with best fit line.  Designated by test temperature.  X and Y labels pulled from \code{data} header for time and stress respectively.
 }
 \item{TestTemp}{
    A vector of the creep test temperatures (from \code{data})
 }
 \item{CreepModelParam}{
    A vector of the creep model parameters for each creep test temperature.  This output only occurs when \code{Tref} is not defined as input.
 }
 \item{RefCreepModelParam}{
    The reference creep model parameter at a test temperature of interest \code{Tref} from which rupture times at other creep test temperature levels are evaluated.  This output only occurs when \code{Tref} is defined as input.
 }
 \item{RuptureTimeatRefCreepModelParam}{
    The rupture time vector of all creep test temperature levels when reference creep model parameter at a test temperature of interest \code{Tref} is defined.  This output only occurs when \code{Tref} is defined as input.
 }
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Reuel Smith,  (2022-2023)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Consider the creep rupture data from a series of tests on Cr-Mo-V steel shown in Table 2.13 of Modarres, Amiri, & Jackson 2017.  Find the Larson-Miller, Manson-Haferd, and Sherby-Dorn parameters for test parameters 900 degrees F at 70 ksi.  Then find the rupture time if the temperature is held at 1100 degrees F at 70 ksi.
creep.analysis(datCreep1,1,list(C = 22, Strace = 70, Tref = 900),2)
creep.analysis(datCreep1,2,list(Strace = 70, Tref = 900,temp_a = 100.13, log10t_a = 18),2)
creep.analysis(datCreep1,3,list(Strace = 70, Tref = 900,Q = 460000),2)


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
