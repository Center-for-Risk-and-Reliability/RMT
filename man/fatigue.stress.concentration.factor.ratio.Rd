\name{fatigue.stress.concentration.factor.ratio}
\alias{fatigue.stress.concentration.factor.ratio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fatigue Stress Concentration Factor Ratio Calculator
}
\description{
This calculator approximates the fatigue notch factor ratio based on the ultimate tensile strength.
}
\usage{
fatigue.stress.concentration.factor.ratio(Su, material, stressunits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Su}{
  Ultimate tensile strength, \eqn{S_u} in ksi or MPa.
}
  \item{material}{
  The name of the material.  Current options for this calculator are limited to "Steel", "Aluminum", and "Magnesium".
}
  \item{stressunits}{
  Defines the units for the stress amplitude: \code{1} for \emph{'MPa'} and \code{2} for \emph{'ksi'}.  \emph{'MPa'} will be used as default if this heading is left empty.
}
}
\details{
The fatigue notch factor ratio is expressed as,
\deqn{(K_f ' - 1)/(K_f - 1)}
where \eqn{K_f} is the fatigue notch factor and \eqn{K_f '} the fatigue notch factor at 1000 cycles.  This ratio is evaluated based on data presented in Bannantine [1] Figure 4.10.  The assumed fit for this data is a fourth degree polynomial with a zero intercept, however because the ratio is bounded between 0 and 1 (where \eqn{(K_f ' - 1)/(K_f - 1) = 0} when \eqn{K_f '=1} and \eqn{(K_f ' - 1)/(K_f - 1) = 1} when \eqn{K_f '= K_f}) any ratio greater than 1 from this polynomial fit resets at 1.
}
\value{
Calculator only returns the fatigue notch factor ratio \eqn{(K_f ' - 1)/(K_f - 1)}.
}
\references{
Banantine, J., Corner, J., & Handrock, J. (1990). Fundamentals of Metal Fatigue Analysis. Prentice Hall
}
\author{
Reuel Smith (2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1: For steel and a given ultimate strength of 114 ksi:
fatigue.stress.concentration.factor.ratio(114, "Steel", 2)

# Example 2: For Aluminum alloy and a given ultimate strength of 600 MPa:
fatigue.stress.concentration.factor.ratio(600, "Aluminum", 1)
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
