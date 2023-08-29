\name{lifestress.PHWbl.LSQest}
\alias{lifestress.PHWbl.LSQest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-Squares PH-Weibull Life-Stress Estimator
}
\description{
Computes least-squares estimate of the Proportional-Hazard Weibull Life-Stress model parameters.
}
\usage{
lifestress.PHWbl.LSQest(pp, interact_stress)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pp}{
The full Weibull probability plot output list from the \code{\link{probplot.wbl}} function.
}
  \item{interact_stress}{
A status which determines whether or not  the estimator will consider interactions between the various stresses or not.  If evaluation checks for interactions \code{interact_stress}, then set to 1, otherwise set to 0 (default).
}
}
\details{
The least-squares estimate of the proportional hazard Weibull life-stress distribution model is based primarily of the stress-based probability plot of the entered data (usually in the \code{probplot} subfunction).  The life-stress parameters are estimated by the relationship linear analysis between stress and the parameter that is based on life (that is \eqn{\alpha}.

The proportional hazard Weibull distribution is positively bounded (\eqn{x \geq 0}) consisting of a Weibull shape parameter \eqn{\beta}, and \eqn{S} and \eqn{\theta} are stress and parameters (respectively) of a modifier funciton \eqn{g(S;\theta)}.

  PDF:
  \deqn{f(x) = \beta x^{\beta - 1} \exp[\Sigma_{j=0}^m \theta_j S_j - x_j^\beta \exp(\Sigma_{j=0}^m \theta_j S_j)]}

  CDF:

  \deqn{F(x) = 1 - \exp[- x_j^\beta \exp(\Sigma_{j=0}^m \theta_j S_j))]}

  Modifier Function:

  \deqn{g(S;\theta) = \exp(\Sigma_{j=0}^m \theta_j S_j)}

The proportional hazard function is formed by multiplying the baseline hazard function (for a typical Weibull distribution) by the modifier funciton, from which we get the relations for PDF and CDF.
}
\value{
Outputs a printed summary of the life-stress distribution model, its least-squares estimate, and average life at each stress level.  The list output includes the following:
\item{S}{
  The vector or matrix of accelerated stresses
}
\item{L}{
  The average life at each stress level
}
\item{LSQ}{
  The least squares estimate for the model parameters (life distribution parameter leads followed by the life-stress parameters
)
}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2022)
}
\note{
This function \strong{only} applies for the proportional hazard Weibull model.  Proportional hazard models with other distributions will be available in the next version.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Parametric PH-Weibull model with two continuous covariates (inverse temperature K^-1 and inverse relative humidity)
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
invtemp <- c(0.002544529, 0.002544529, 0.002544529, 0.002544529, 0.002544529, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861, 0.002832861)
invRH <- c(1.666666667, 1.666666667, 1.666666667, 1.666666667, 1.666666667, 1.176470588, 1.176470588, 1.176470588, 1.176470588, 1.176470588, 1.666666667, 1.666666667, 1.666666667, 1.666666667, 1.666666667)
mat1 <- matrix(c(time, cens, invtemp,invRH), nrow = length(time), ncol = 4, byrow = FALSE,dimnames = list(NULL,c("time (hrs)","cens","invtemp","invrelhumidity")))
lifestress.PHWbl.LSQest(probplot.wbl(mat1,"Blom","hours")[[1]], 0)

## Example 2: Parametric PH-Weibull model with two continuous covariates (inverse temperature K^-1 and inverse relative humidity), but in this instance, consider the interactions between the two stresses
lifestress.PHWbl.LSQest(probplot.wbl(mat1,"Blom","hours")[[1]], 1)

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
