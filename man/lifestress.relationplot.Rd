\name{lifestress.relationplot}
\alias{lifestress.relationplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Life-Stress Relationship Plot Generator
}
\description{
Post-processes the relationship plot for a given stress-life model and output parameters (least-squares and maximum-likelihood estimate).
}
\usage{
lifestress.relationplot(data, ls, dist, pp, Smin, Smax, SUse, therm, confid, Llab, Slab)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 consists of the failure (or primary event) data and the right censored data in a given set.  Data must consist of real and finite numbers.
  }
  \item{
  Column 2 gives the censored status of the data in column 1.  Demark failed data as \eqn{1} and right-censored data as \eqn{0}.
  }
  \item{
  Column(s) 3 (and on for multi-type stress) states the stress values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).  Column 3 is not necessary if all data is from the same stress level however.
  }
}
}
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
  \item{Smin}{
The minimum value for stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).
}
  \item{Smax}{
The maximum value for stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).
}
  \item{SUse}{
The use-level stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).
}
  \item{therm}{
Designation of the primary thermal stress in the original data.  If it is the first column after the censored column then the value is 1 and the alternate is 2 (this is set to default).  If it is the second column after the censored column then the value is 2 and the alternate is 1.  Only necessary to input for dual-stress life-stress models.
}
  \item{confid}{
The confidence in the maximum likelihood parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{Llab}{
Written label used for the life axis.
}
  \item{Slab}{
Written label used for the charateristic stress axis or axes.  When the Arrhenius or Eyring life-stress model is selected, this label defaults to \eqn{K^{-1}} where \eqn{K} is the temperature notation for Kelvin.  Dual-stresses must be entered as vectors (Example: For "Temperature" and "Humidity" enter \code{c("Temperature", "Humidity")}).
}
}
\details{
The characteristic relationship plot compares a set of given stresses to the characteristic life defined by the life-stress model \code{ls}.  The point of reference for the relationship plot is the use stress.
}
\value{
The function provides the relationship plot for the given input in addition to the output list of life values at use-stress for the least-squares estimate (output [[1]]) and the maximum likelihood estimate (output [[2]]).
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single Stress value (Temperature in Kelvin)
time <- c(248, 456, 528, 731, 813, 965, 972, 1528, 164, 176, 289, 319, 386, 459, 92, 105, 155, 184, 219, 235)
cens <- rep(1,length(time))
temp <- c(406, 406, 406, 406, 406, 406, 406, 406, 436, 436, 436, 436, 436, 436, 466, 466, 466, 466, 466, 466)
mat1 <- matrix(c(time, cens, temp), nrow = length(time), ncol = 3, byrow = FALSE,dimnames = list(NULL,c("time (hrs)","cens","temp")))
lifestress.relationplot(mat1,"Arrhenius","Weibull",probplot.wbl(mat1,"Blom","Time (hours)")[[1]],200,490,293.15,0,0.95,"Time (hours)","Degrees (Kelvin)")

## Example 2: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat2 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE,dimnames = list(NULL,c("time (hrs)","cens","temp","Relative humidity")))
lifestress.relationplot(mat2,"TempHumidity","Weibull",probplot.wbl(mat2,"Blom","Time (hours)")[[1]],c(250,0.5),c(400,.9),c(293.15,0.2),0,0.95,"Time (hours)",c("Degrees (Kelvin)","Relative Humidity"))
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
