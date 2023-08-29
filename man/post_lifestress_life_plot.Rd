\name{post_lifestress_life_plot}
\alias{post_lifestress_life_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior Life-Stress Relationship Plot
}
\description{
Post-processes the relationship plot for a given stress-life model and output parameters (Bayesian parameter estimate).
}
\usage{
post_lifestress_life_plot(post_params, ls, confid, Xlab, Ylab, Smin, Smax, S0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{post_params}{
  The mean and standard deviation of the posterior parameter markov chains as well as the log-likelihood (Output[[3]] from \code{lifestress.BAYESest}).
}
  \item{ls}{
The life-stress model definition for analysis.  Choose from the list of life-models outlined in \code{\link{lifestress.select}}.
}
  \item{confid}{
The confidence in the Bayesian parameter estimates entered as a probability between \eqn{0} and \eqn{1}.  Default is \eqn{0.95} (or \eqn{95\%})
}
  \item{Xlab}{
Written label used for the charateristic stress axis or axes.  When the Arrhenius or Eyring life-stress model is selected, this label defaults to \eqn{K^-1} where \eqn{K} is the temperature notation for Kelvin.  Dual-stresses must be entered as vectors (Example: For "Temperature" and "Humidity" enter \code{c("Temperature", "Humidity")}).
}
  \item{Ylab}{
Written label used for the life axis.
}
  \item{Smin}{
The minimum value for stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).
}
  \item{Smax}{
The maximum value for stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).
}
  \item{S0}{
The use-level stress to be used for the relationship plot (if dual stress enter as a vector \code{c(Stress1,Stress2)}).  This is an optional input but if left out, then the output values will vary (see below).
}
}
\details{
The characteristic relationship plot compares a set of given stresses to the characteristic life defined by the life-stress model \code{ls}.  The point of reference for the relationship plot is the use stress.
}
\value{
The function provides the relationship plot for the given input in addition to the following depending on if use stress \code{S0} is defined.  If \code{S0} is defined:

  \item{plotout2}{
  The relationship plot consisting of the upper and lower confidence band of the
  }
  \item{comp2 }{Description of 'comp2'}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Reuel Smith (2020-2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


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
