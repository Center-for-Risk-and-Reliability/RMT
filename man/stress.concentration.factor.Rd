\name{stress.concentration.factor}
\alias{stress.concentration.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stress Concentration Factor Calculator
}
\description{
Computes the stress concentration factor \eqn{K_t} for a specified geometry and set of dimensions.
}
\usage{
stress.concentration.factor(dimensions, geometry, loadtype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dimensions}{
Standard dimensions of the sample geometry based on several of the following: Width \code{W}, thickness \code{t}, semicircular notch radius \code{r}, U-notch depth \code{d}, and/or long diameter \code{d}.  These dimensions must be entered as a list.
}
  \item{geometry}{
Definition of the sample geometry from which a stress concentration factor is sought from.  The current list of geometries is as follows:

\strong{Rectangular Bar or Plate options}
  \itemize{
      \item \strong{RB1. Rectangular Bar with a Semi-circle Edge Notch} - Enter as \code{"rect_1semicirc_edge"} with dimensions width \code{W} and radius \code{r}
      \item \strong{RB2. Rectangular Bar with Opposite Semi-circle Edge Notches} - Enter as \code{"rect_2semicirc_edge"} with dimensions width \code{W} and radius \code{r}
      \item \strong{RB3. Rectangular Bar with a U-Notch} - Enter as \code{"rect_1U_edge"} with dimensions width \code{W}, radius \code{r}, and U-depth \code{d}
      \item \strong{RB4. Rectangular Bar with Opposite Edge U-Notches} - Enter as \code{"rect_2U_edge"} with dimensions width \code{W}, radius \code{r}, and U-depth \code{d}
      \item \strong{RB7. Plate with a Circular Symmetric Hole} - Enter as \code{"rect_circ_symm_hole"} with dimensions width \code{W}, short diameter \code{a}, and long diameter \code{d}
  }
}
  \item{loadtype}{
The description of the load on the sample: Enter \code{1} for tension (default), enter \code{2} for bending, and enter \code{3} for traverse bending.  (See note on options 2 and 3)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Computes the stress concentration factor \eqn{K_t}.
}
\references{
Pikley, W., Pikley, D., & Bi, Z. (2020). Peterson's Stress Concentration Factors (4th Edition). Wiley
}
\author{
Reuel Smith (2022)
}
\note{
The full list of geometries is currently incomplete.  Will add more options in later updates.  Options for \code{loadtype} are currently limited to tension.  The other two will be added in a later update.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
stress.concentration.factor(list(W=100,r=10),"rect_2semicirc_edge")

stress.concentration.factor(list(W=100,r=10),"rect_1semicirc_edge")

stress.concentration.factor(list(r = 20, W = 100),"rect_2semicirc_edge")

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
