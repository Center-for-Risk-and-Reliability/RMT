\name{probplot.logist}
\alias{probplot.logist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic Probability Plot
}
\description{
Plots the Logistic probability plot for singular and multiple stress levels based on an input data set (a matrix) and a given plotting position.
}
\usage{
probplot.logist(data,pp,xlabel1="X")
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
  \item{pp}{
Name of the plotting position to select for computation.  See \code{\link{plotposit.select}} for plotting position definitions.
}
  \item{xlabel1}{
Label for the x-axis for the plot.  If left blank, the default label is "X".
}
}
\details{
The Logistic distribution is unbounded (\eqn{-\infty < x < \infty}) and consists of a logistic location parameter \eqn{\mu} and a logistic scale parameter \eqn{\sigma}.

  PDF:
  \deqn{f(x) = \frac{\exp(\frac{x - \mu}{\sigma})}{\sigma [1+\exp(\frac{x - \mu}{\sigma})]^2}}

  CDF:

  \deqn{F(x) = \frac{\exp(\frac{x - \mu}{\sigma})}{1+\exp(\frac{x - \mu}{\sigma})}}
}
\value{
The output, in addition to the logistic probability plot, includes the following for each stress group outlined in the \code{data} input:

\item{Stress Group}{
Lists the stress level based on the grouped data.  These will always be the first of each group of three on the output list (i.e. [[1]], [[4]], [[7]], etc.).
  }
  \item{logistresults}{
  The least squares estimates for the logistic location parameter \eqn{\mu} and a logistic scale parameter \eqn{\sigma}.  These will always be the second of each group of three on the output list (i.e. [[2]], [[5]], [[8]], etc.).
}
\item{R2}{
  Coefficient of determination \eqn{R^2}.  These will always be the third of each group of three on the output list (i.e. [[3]], [[6]], [[9]], etc.).
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2023)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single data set (no stress distinction).  NOTE: Data must include a stress column because the group stress definition is one of the output.
mat0 <- matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8)
temp <- rep(373,length(mat0[,1]))
mat1 <- matrix(c(mat0[,1], mat0[,2], temp), nrow = length(mat0[,1]), ncol = 3, byrow = FALSE)
probplot.logist(mat1, "Blom", "Time (hours)")

## Example 2: Single Stress value (Temperature in Kelvin)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
mat2 <- matrix(c(time, cens, temp), nrow = length(time), ncol = 3, byrow = FALSE)
probplot.logist(mat2, "Blom", "Time (hours)")

## Example 3: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat3 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
probplot.logist(mat3, "KaplanMeier", "Time (hours)")
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
