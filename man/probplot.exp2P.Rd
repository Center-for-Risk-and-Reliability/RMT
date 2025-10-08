\name{probplot.exp2P}
\alias{probplot.exp2P}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Parameter Exponential Probability Plot
}
\description{
Plots the two-parameter exponential probability plot for singular and multiple stress levels based on an input data set and a given plotting position.  The default plot is based on least squared (LSQ) estimates for the parameters.
}
\usage{
probplot.exp2P(data, pp, xlabel1, confid, nobounds, MLE_i, CDFrangesetting, stressunit1, stressunit2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The sample set of event or time data to be plotted.  The standard format is a matrix formed as follows:
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
Data may now be simply entered as a vector where the assumption is that all entered data is considered failed and irrespective of any stress value.
  }
  \item{pp}{
Name of the plotting position to select for computation.  See \code{\link{plotposit.select}} for plotting position definitions.  The default is \code{"Blom"}
  }
  \item{xlabel1}{
Label for the x-axis for the plot.  If left blank, the default label is "X".
  }
  \item{confid}{
The confidence in the best fit probability plot curve entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
  }
  \item{nobounds}{
(NEW) Option allows you to turn off the confidence interval plotting where the default is \code{NULL}.  Enter \code{1} to deactivate the confidence intervals.
  }
  \item{MLE_i}{
(NEW) Activates calculation of the maximum likelihood estimates (MLE) for the parameters and plot for the probability plot.  Default is \code{NULL} for LSQ.  Enter \code{1} to compute the MLE.
  }
  \item{CDFrangesetting}{
(NEW) Sets the failure probability range for the y-axis.  Default is \code{1} for 1\% to 99\%.  Enter  \code{2} for 0.1\% to 99.9\%.
  }
  \item{stressunit1}{
(NEW) Label for the first stress unit.  If left blank, the default label is "unit".
  }
  \item{stressunit2}{
(NEW) Label for the second stress unit (if applicable).  If left blank, the default label is "unit".
  }
}
\details{
  Two-parameter exponential distribution is positively bounded (\eqn{x \geq 0}) consisting of a two-parameter exponential rate parameter \eqn{\sigma} and a two-parameter exponential location parameter \eqn{\theta}.

  PDF:
  \deqn{f(x) = \frac{1}{\sigma} \exp[- \frac{x - \theta}{\sigma}]}

  CDF:

  \deqn{F(x) = 1 - \exp[- \frac{x - \theta}{\sigma}]}
}
\value{
Output consists of a table of parameter estimates by stress as well as the confidence interval for each parameter if \code{MLE_i=1}.  The remaining output list consists of the following:
  \item{output}{
Lists the primary output of estimates by stress.  Each group includes the following by default (LSQ):
\itemize{
  \item{
  Stress Level - The unitless value or values of stress
  }
  \item{
  Parameter Estimates - The LSQ estimates for the two-parameter exponential rate parameter \eqn{\sigma} and the two-parameter exponential location parameter \eqn{\theta}.
  }
  \item{
  SSE - The sum of square error (based on the linearlized CDF.
  }
  \item{
  Coefficient of Determination (R2) - The  Coefficient of determination \eqn{R^2}.
  }
}
If MLE, the grouping will be as follows:
\itemize{
  \item{
  Stress Level - The unitless value or values of stress
  }
  \item{
  Parameter Estimates - The MLE estimates for the two-parameter exponential rate parameter \eqn{\sigma} and the two-parameter exponential location parameter \eqn{\theta}
  }
  \item{
  CI - The confidence intervals for each MLE parameter estimate (based on \code{confid}).
  }
  \item{
  SSE - The sum of square error (based on the linearlized CDF.
  }
  \item{
  loglikelihood - The  loglikelihood value based on the MLE point estimate.
  }
  \item{
  likelihood - The likelihood value based on the MLE point estimate.
  }
  \item{
  AIC - The Akaike Information Criterion value.
  }
  \item{
  BIC - The Bayesian Information Criterion value.
  }
}
  }
    \item{summary.nonparametric}{
Lists the nonparametric output table for each stress value:Column 1 \emph{failure event data or time}, is  Column 2 is the \emph{failure probability}, Column 3 is the \emph{reliability}, Column 4 is the \emph{hazard rate}, Column 5 is the \emph{cumulative hazard}, and Column 6 is the \emph{probability density}.
  }
    \item{prob_plot}{
Generated probability plot.
  }
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2025)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single data set (no stress distinction).
mat0 <- matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8)
temp <- rep(373,length(mat0[,1]))
mat1 <- cbind(mat0,temp)
probplot.exp2P(mat1, "Blom", "Time (hours)",stressunit1 = "K")

## Example 2: Single Stress value (Temperature in Kelvin)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
mat2 <- cbind(time, cens, temp)
probplot.exp2P(mat2, "Blom", "Time (hours)",stressunit1 = "K")

## Example 3: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat3 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
probplot.exp2P(mat3, "Blom", "Time (hours)",stressunit1 = "K", stressunit2 = "RH")

# Example 4: Repeat Example 3 but with the MLE estimates.
probplot.exp2P(mat3, pp = "Blom", xlabel1 = "Time (hours)",MLE_i = 1,stressunit1 = "K", stressunit2 = "RH")

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
