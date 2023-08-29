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
Full line of life-stress distribution parameters.  The output from any of the life-stress esimator functions (\code{lifestress.LSQest}, \code{lifestress.MLEest}, or \code{lifestress.BAYESest}).
}
  \item{lsm}{
The life-stress model definition for analysis.  Choose from the list of life-models outlined in \code{\link{lifestress.select}}.
}
  \item{Sacc}{
The accelerated stress for the problem.  If dual or multiple stress, enter as a vector starting with temperature based stress.
}
  \item{Suse}{
The use stress for the problem.  If dual or multiple stress, enter as a vector starting with temperature based stress.
}
}
\details{
Acceleration factor \eqn{AF} calculations are based on the given stress-life equations \eqn{L(S)} and the standard notation where,

\deqn{AF = L(S_use)/L(S_acc)}
}
\value{
Returns the acceleration factor for this stress scenario.
}
\references{
%% ~put references to the literature/web site here ~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (params, lsm, Sacc, Suse)
{
    lsfcn <- lifestress.select(lsm)[[1]]
    if (length(Suse) == 1 | lsm == "MultiStress") {
        AF <- lsfcn(params, Suse)/lsfcn(params, Sacc)
    }
    if (length(Suse) == 2 & !(lsm == "MultiStress")) {
        AF <- lsfcn(params, Suse[1], Suse[2])/lsfcn(params, Sacc[1],
            Sacc[2])
    }
    return(AF)
  }
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
