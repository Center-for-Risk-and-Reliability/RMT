\name{corrosion.pitting}
\alias{corrosion.pitting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pitting Corrosion Calculator
}
\description{
Analyzes pit growth data due to pitting corrosion using various chemical properties of the material and assuming the Kondo-Wei model to define pit size.
}
\usage{
corrosion.pitting(data, r_0, r_cr, matproperties, units)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Tabular pit corrosion data with labeled headers.  Column 1 consists of the time data, Column 2 is the pit radius or growth data, Column 3 is the pit name ("Pit 1", "Pit 2", etc.), Column 4 is the temperature in Kelvin (may be in Rankine, but calculator will convert units to Kelvin), and optional Column 5 is the pH value of the pit.
  }
  \item{r_0}{
  Initial pit radius (0 is default)
  }
  \item{r_cr}{
  Final or critical pit radius.
  }
  \item{matproperties}{
A list made up of the properties of the material under study: Molecular weight \eqn{M} (\code{M}) in g/mol, material density \eqn{\rho} (\code{rho}) in \eqn{\text{g}/\text{m}^3}, and number of valence electrons \eqn{n} (\code{n}).
}
  \item{units}{
Defines the units for the measurements: \code{1} for Metric or SI (temperature \emph{'Kelvin'}) and \code{2} for English (temperature \emph{'Rankine'}).  Metric will be used as default if this heading is left empty.
}
}
\details{
The pitting corrosion calculator performs its analysis based on the Kondo-Wei model that defines pit size as a hemispherical volume of radius \eqn{r}.  The volumetric pitting corrosion rate \eqn{i_{corr}} is defined under this model as,

\deqn{i_{corr}=2\pi r^2\frac{dr}{dt}=\frac{M I_{p_{0}}}{nF\rho}\exp\left(-\frac{E_a}{k_B T}\right)}

where \eqn{F} is the Faraday constant \eqn{96,514} C/mol, \eqn{k_B} is the Boltzmann constant \eqn{8.617\times 10^{-5}} eV/K, \eqn{E_a} is the activation energy, \eqn{T} is the temperature in Kelvin, and \eqn{I_{p_{0}}} is the pitting current coefficient in C/s defined by \eqn{I_{p_{0}} = \left(1+p\right)\times 6.5\times10^{-5}} where \eqn{p} is the number of particles in a group. We can further define time of corrosion \eqn{t} from an initial pit size of radius \eqn{r_0} to a final or critical pit size of radius \eqn{r_{f}} as,

\deqn{t=\frac{2\pi}{3}\frac{nF\rho}{M I_{p_{0}}}\left(r_{f}^3 - r_{0}^3\right)\exp\left(\frac{E_a}{k_B T}\right)}

The pitting corrosion rate may also be defined in terms of temperature and pH where,

\deqn{\frac{2 \pi}{3}\frac{nF\rho}{M I_{p_{0}}}=A\times pH}

where \eqn{A} is a model parameter in time per volume.  We can redefine corrosion rate \eqn{i_{corr}} and corrosion time \eqn{t} respectively under this relation as,

\deqn{i_{corr}=2\pi r^2\frac{dr}{dt}=\frac{2\pi}{3A\times pH}\exp\left(-\frac{E_a}{k_B T}\right)}

\deqn{t=A\times pH \times\left(r_{f}^3 - r_{0}^3\right)\exp\left(\frac{E_a}{k_B T}\right)}

}
\value{
The pitting corrosion calculator will produce the following output based primarily on the data provided.
    \item{corroutputtable}{
        A pit-by-pit output table consisting of the Kondo-Wei model parameters, pitting corrosion time from \eqn{r_0} to \eqn{r_{f}}, the pitting corrosion rate, and the coefficient of determination \eqn{R^2}.  The model parameters produced will vary based on the input.  If \emph{temperature} and \emph{material parameters} are provided then \eqn{I_{p_{0}}} and \eqn{E_a} will be calculated.  If \emph{temperature} and \emph{pH} are provided, then \eqn{A} and \eqn{E_a} will be calculated.  If all three \emph{temperature}, \emph{pH}, and \emph{material parameters} are provided then \eqn{A}, \eqn{I_{p_{0}}}, and \eqn{E_a} will be calculated.
    }
    \item{stresstimetable}{
        This is a standard time-to-failure table with accompanying stress values.  It may be used for ALT analysis.
    }
    \item{corrplot}{
        Plot of the test time vs. pit radius growth for all pits from the input data.
    }
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability

Y. Kondo, Y. (1989) "Prediction of Fatigue Crack Initiation Life Based on Pit Growth," Corrosion, Vol. 45, No. 1, pp. 7-11.

Kondo, Y.  & Wei, R.P. (1989) "Approach on Quantitative Evaluation of Corrosion Fatigue Crack Initiation Condition," International Conference on Evaluation of Materials Performance in Severe Environments, EVALMAT 89, Vol. 1, Iron and Steel Institute of Japan, pp. 135-142.


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
## Example 1: Download the following data

datPitexample  <- read.csv("https://raw.githubusercontent.com/Center-for-Risk-and-Reliability/RMT/main/CSVExampleData/Pitdataexample.csv")

## This is pit growth data as related to operation temperature and known pH value.  Calculate the Kondo-Wei parameters, corrosion time, and corrosion rate for this data where the critical radius is given as 0.05 mm.

corrosion.pitting(datPitexample, r_0 = 0, r_cr = 0.05, units = 1)

## Assume that we now know this to be from an Aluminum alloy material where density is 2.7e6 g/m^3, molecular weight is 27 g/mol, and valence number is 3.  Now compute the Kondo-Wei parameters, corrosion time, and corrosion rate based on this information.

corrosion.pitting(datPitexample, r_0 = 0, r_cr = 0.05, matproperties = list(rho = 2.7e6, M = 27, n = 3),units = 1)
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
