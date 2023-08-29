\name{hysteresisloop.plot}
\alias{hysteresisloop.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hysteresis Loop Plot
}
\description{
Generates a plot of the strain-life hysteresis loop for loading and load cycling.
}
\usage{
hysteresisloop.plot(E,K,n,stressunits,loadconditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{
      Modulus of elasticity \eqn{E} of the material (in MPa or ksi).  This needs to be in the same units as that of stresses and ranges in \code{loadconditions}.
  }
  \item{K}{
      Cyclic strength coefficient \eqn{K'}.
  }
  \item{n}{
      Cyclic strain hardening component \eqn{n'}.
  }
  \item{stressunits}{
      Defines the units for the stress amplitude: \code{1} for \emph{'MPa'} and \code{2} for \emph{'ksi'}.  \emph{'MPa'} will be used as default if this heading is left empty.
  }
  \item{loadconditions}{
      A list containing one of six possible inputs:
      \itemize{
          \item \strong{Maximum or minimum stress}
          (Enter as \code{maxstress} or \code{minstress} respectively).  Enter one of these when loading is fully reversed (\eqn{R=-1}).
          \item \strong{Stress Range}
          (Enter as \code{stressrange}).  When loading is not fully reversed enter the stress loading range as a vector \code{c(maxtress,minstress)}.
          \item \strong{Maximum or minimum strain}
          (Enter as \code{maxstrain} or \code{minstrain} respectively).  Enter one of these when loading is fully reversed (\eqn{R=-1}).
          \item \strong{Strain Range}
          (Enter as \code{strainrange}).  When loading is not fully reversed enter the strain loading range as a vector \code{c(maxstrain,minstrain)}.
      }
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list consisting of the following:
  \item{hysteresisloop}{
      The plot of the hysteresis loop
  }
  \item{maxstrain}{
  The maximum strain as given or as calculated in the case of a stress-based input.
  }
  \item{minstrain}{
  The minimum strain as given or as calculated in the case of a stress-based input.
  }
  \item{maxstress}{
  The maximum stress as given or as calculated in the case of a strain-based input.
  }
  \item{minstress}{
  The minimum stress as given or as calculated in the case of a strain-based input.
  }
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Reuel Smith, (2022)
}
\note{
The input \code{loadconditions} is left in list form as another input will be made available in a future update.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1: A material has properties stated in Example 2.3 E = 210000 MPa, n = 0.2, and K = 1200 MPa.  Assuming a fully reversed cycle where delta strain is 0.04 we can enter either maxstrain = 0.02 or minstrain = -0.02 as our load condition.
hysteresisloop.plot(210000,1200,0.2,1,list(maxstrain = 0.02))

# Example 2: Similarly we can enter a max stress as 534 MPa or min stress of -534 MPa for the same example.
hysteresisloop.plot(210000,1200,0.2,1,list(maxstress = 534))

# Example 3: If instead we set a load from -500 to 540 MPa.
hysteresisloop.plot(210000,1200,0.2,1,list(stressrange = c(-500,540)))
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
