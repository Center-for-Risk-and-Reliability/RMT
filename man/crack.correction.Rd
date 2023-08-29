\name{crack.correction}
\alias{crack.correction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crack Correction Factor Calculator
}
\description{
Computes the crack correction factor (\eqn{f(g)} or \eqn{Y}) for a specified geometry and set of dimensions.
}
\usage{
crack.correction(dimensions, crackgeometry)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dimensions}{
        A list of standard dimensions of the sample geometry based on several of the following: Width \code{W}, thickness \code{t}, semicircular notch radius or depth (or initial center crack size \eqn{2 r}) \code{r}, and/or initial crack size beyond \code{r}, \code{a_0}.  These dimensions must be entered as a list.
  }
  \item{crackgeometry}{
       Definition of the sample geometry from which a crack correction factor is sought from if one is not provided in \code{data}.  The current available geometries are as follows:
       \itemize{
      \item \strong{Center cracked plate} - Enter as \code{"center_2a"} with dimensions width \code{W} and half-center crack or notch length  \code{r}
      \item \strong{Single Edge cracked plate} - Enter as \code{"edge_single_1a"} with dimensions width \code{W} and radius or notch length \code{r}
      \item \strong{Double  Edge cracked plate} - Enter as \code{"edge_double_2a"} with dimensions width \code{W} and radius or notch length \code{r}
  }
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The crack correction factor (\eqn{f(g)} or \eqn{Y}).
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability

Banantine, J., Corner, J., & Handrock, J. (1990). Fundamentals of Metal Fatigue Analysis. Prentice Hall
}
\author{
Reuel Smith (2022)
}
\note{
The full list of geometries is currently incomplete.  Will add more options in later updates.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1
crack.correction(list(r = 0.5, W = 100, t = 5),"edge_single_1a")

# Example 2
crack.correction(list(r = 0.5, W = 100, t = 5),"edge_double_2a")

# Example 3
crack.correction(list(r = 0.5, W = 100, t = 5),"center_2a")
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
