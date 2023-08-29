\name{notch.SN.diagram}
\alias{notch.SN.diagram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stress-Life Diagram with Notch Effect
}
\description{
Generates the stress-life (or S-N) diagram based on various stress and fatigue life input data as well as notch effect.
}
\usage{
notch.SN.diagram(dimensions, geometry, stressunits, loadconditions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dimensions}{
Standard dimensions of the sample geometry based on several of the following: Width \code{W}, thickness \code{t}, semicircular notch radius \code{r}, U-notch depth \code{d}, and/or long diameter \code{d}.  These dimensions must be entered as a list.
}
  \item{geometry}{
Definition of the sample geometry from which a stress concentration factor is sought from.  The current available geometries are listed in \code{\link{stress.concentration.factor}}.
}
  \item{stressunits}{
Defines the units for the stress amplitude: \code{1} for \emph{'MPa'} and \code{2} for \emph{'ksi'}.  \emph{'MPa'} will be used as default if this heading is left empty.
}
  \item{loadconditions}{
A list of load conditions and/or material properties that may be required to compute notch effects.  Includes: (1) \code{Su} the ultimate strength \eqn{S_u}, (2) \code{Se} the fatigue limit \eqn{S_e}, (3) \code{Sa} the stress amplitude \eqn{S_a}, and/or (4) \code{a} the material constant needed to compute the notch sensitivity factor \eqn{q} in either inches or mm.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Produces an S-N figure made up of two curves: (1) the unnotched curve and (2) the notched curve.  Also produces a list of the following output:
  \item{Kt}{
      The stress concentration factor \eqn{K_t} evaluated by \code{\link{stress.concentration.factor}}.
  }
  \item{q}{
      The notch sensitivity factor \eqn{q}.
  }
  \item{Kf}{
      The fatigue stress concentration factor \eqn{K_f}.
  }
  \item{Kfp}{
      The fatigue notch factor at 1000 cycles \eqn{K_f '}.  This is evaluated based on \eqn{K_f} and the the fatigue notch factor ratio \eqn{(K_f ' - 1)/(K_f - 1)} found by \code{\link{fatigue.stress.concentration.factor.ratio}}.
  }
  \item{unnotched_S}{
      The list of alternating stress values at \eqn{1}, \eqn{1,000}, and \eqn{1,000,000} fatigue cycles for an unnotched sample.
  }
  \item{notched_S}{
      The list of alternating stress values at \eqn{1}, \eqn{1,000}, and \eqn{1,000,000} fatigue cycles for a notched sample.
  }
  \item{Snet}{
      The net stress \eqn{S_{net}} of a sample based on its geometry.
  }
  \item{A}{
    The estimate for stress-life parameter \eqn{A} based on notch effect.
  }
  \item{b}{
    The estimate for stress-life parameter \eqn{b} based on notch effect.
  }
  \item{Nf}{
    The estimate notched stress-life at \eqn{S_{net}}.
  }
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability

Banantine, J., Corner, J., & Handrock, J. (1990). Fundamentals of Metal Fatigue Analysis. Prentice Hall
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
# Example 1: Example 2.7 from the textbook
notch.SN.diagram(list(W = 100, r = 10, t = 5),"rect_2semicirc_edge",1,list(Su = 790, Kfratio = 0.3, Sa = 150000))

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
