\name{lifestress.LSQest}
\alias{lifestress.LSQest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-Squares Life-Stress Estimator
}
\description{
Computes least-squares estimate of Life-Stress model parameters for a given set of data assigned a specific life distribution.
}
\usage{
lifestress.LSQest(ls, dist, pp, therm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ls}{
The life-stress model definition for analysis.  Choose from the list of life-models outlined in \code{\link{lifestress.select}}.
}
  \item{dist}{
  The distribution that represents time or life.
  \enumerate{
  \item \strong{Weibull} (Enter as "Weibull".  See more on \code{\link{probplot.wbl}}.)
  \item \strong{Three Parameter Weibull} (Enter as "3PWeibull".  See more on \code{\link{probplot.wbl3P}}.)
  \item \strong{Exponential} (Enter as "Exponential".  See more on \code{\link{probplot.exp}}.)
  \item \strong{Two Parameter Exponential} (Enter as "2PExponential".  See more on \code{\link{probplot.exp2P}}.)
  \item \strong{Normal} (Enter as "Normal".  See more on \code{\link{probplot.nor}}.)
  \item \strong{Lognormal} (Enter as "Lognormal".  See more on \code{\link{probplot.logn}}.)
  }
}
  \item{pp}{
The full probability plot output list from any one of the functions available: (1) \code{\link{probplot.wbl}}, (2) \code{\link{probplot.wbl3P}}, (3) \code{\link{probplot.exp}}, (4) \code{\link{probplot.exp2P}}, (5) \code{\link{probplot.nor}}, and (6) \code{\link{probplot.logn}}.
}
  \item{therm}{
Designation of the primary thermal stress in the original data.  If it is the first column after the censored column then the value is 1 and the alternate is 2 (this is set to default).  If it is the second column after the censored column then the value is 2 and the alternate is 1.  Only necessary to input for dual-stress life-stress models.
}
}
\details{
The least-squares estimate of an entered life-stress distribution model is based primarily of the stress-based probability plot of the entered data (usually in the \code{probplot} subfunction).  The life-stress parameters are estimated by the relationship linear analysis between stress and a parameter that is based on life (for Weibull or three-parameter Weibull that is \eqn{\alpha}, for normal that is \eqn{\mu}, for lognormal that is \eqn{\mu_t}, for exponential that is \eqn{\lambda}, and for two parameter exponential that is \eqn{\sigma}).  The final distribution paramter is obtained by averaging through the multiple stress lines.
}
\value{
Outputs a printed summary of the life-stress distribution model, its least-squares estimate, and average life at each stress level.  The list output includes the following:
%%  If it is a LIST, use
\item{S}{
  The vector or matrix of accelerated stresses
}
\item{L}{
  The average life at each stress lecel
}
\item{LSQ}{
  The least squares estimate for the model parameters (life distribution parameter leads followed by the life-stress parameters)
}
\item{R2}{
  The coefficient of determination.
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
This function can \strong{only} operate for accelerated stress test data.  For step-stress test data use \code{stepstress.LSQest}.  Data matrix \code{mat} may also be loaded from an Excel .CSV file by \code{mat1 <-read.csv("c:/.../Data.csv")}.  This would be used to provide the \code{pp} output needed to run this function (by way of \code{probplot.NAME(mat,"plotposition","Xlabel")}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single Stress value (Temperature in Kelvin)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
mat1 <- matrix(c(time, cens, temp), nrow = length(time), ncol = 3, byrow = FALSE,dimnames = list(NULL,c("time (hrs)","cens","temp")))
lifestress.LSQest("Arrhenius","Weibull",probplot.wbl(mat1,"Blom","hours")[[1]])

## Example 2: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat2 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE,dimnames = list(NULL,c("time (hrs)","cens","temp","Relative humidity")))
lifestress.LSQest("TempHumidity","Weibull",probplot.wbl(mat2,"Blom","hours")[[1]])

## Data(dat) may also be loaded from an Excel .CSV file by dat <-read.csv("c:/.../Data.csv")

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
