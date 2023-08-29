\name{SN.meanstresseffect}
\alias{SN.meanstresseffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stress-Life Mean Stress Effect Calculator
}
\description{
This calculator computes the equivalent stress amplitude \eqn{S_{ar}} for a non-zero mean stress \eqn{S_a} or non-fully reversible scenario.
}
\usage{
SN.meanstresseffect(relationship, Sa, Sm, altvar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relationship}{
The mean stress correction relationship between equivalent stress amplitude \eqn{S_{ar}} for a non-zero mean stress amplitude \eqn{S_a}.
\enumerate{
    \item \strong{Soderberg relationship} - (enter as "Soderberg")
    \deqn{S_a/S_{ar} + S_m/S_y = 1} (\emph{NOTE}: \eqn{S_y} yield stress to be entered in list of \code{altvar})
    \item \strong{Modified Goodman relationship} - (enter as "ModifiedGoodman")
    \deqn{S_a/S_{ar} + S_m/S_u = 1} (\emph{NOTE}: \eqn{S_u} ultimate stress to be entered in list of \code{altvar})
    \item \strong{Morrow relationship} - (enter as "Morrow")
    \deqn{S_a/S_{ar} + S_m/sig_f = 1} (\emph{NOTE}: \eqn{sigf} true fracture strength to be entered in list of \code{altvar})
    \item \strong{Gerber relationship} - (enter as "Gerber")
    \deqn{S_a/S_{ar} + (S_m/S_u)^2 = 1} (\emph{NOTE}: \eqn{S_u} ultimate stress to be entered in list of \code{altvar})
    \item \strong{Smith-Watson-Topper relationship} - (enter as "SWT")
    \deqn{S_a/S_{ar} = ((1 - R)/2)^{0.5}} (\emph{NOTE}: \eqn{R} load ratio to be entered in list of \code{altvar})
    \item \strong{Walker relationship} - (enter as "Walker")
    \deqn{S_a/S_{ar} = ((1 - R)/2)^{(1 - gam)}} (\emph{NOTE}: \eqn{R} load ratio  and \eqn{gam} the fitting constant to be entered in list of \code{altvar})
  }
}
  \item{Sa}{
  The non-zero mean stress amplitude \eqn{S_a}.
}
  \item{Sm}{
  The mean stress \eqn{S_m}.
}
  \item{altvar}{
  Additional input variables for each relationship (see relationship list).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Computes the equivalent stress amplitude \eqn{S_{ar}} for a non-zero mean stress \eqn{S_a}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Reuel Smith, (2022)
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
