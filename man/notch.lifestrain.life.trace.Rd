\name{notch.lifestrain.life.trace}
\alias{notch.lifestrain.life.trace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Life from Life Strain Properties with Notch Effect
}
\description{
Computes the strain-life (in fatigue cycles \eqn{N_f} and reversals) for a given set of strain-life parameters and conditions with consideration for notch effect.
}
\usage{
notch.lifestrain.life.trace(dimensions, geometry, E, b, c, sig_f, eps_f, K, n, stressunits, loadconditions, options)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dimensions}{
Standard dimensions of the sample geometry based on several of the following: Width \code{W}, thickness \code{t}, semicircular notch radius \code{r}, U-notch depth \code{d}, and/or long diameter \code{d}.  These dimensions must be entered as a list.
}
  \item{geometry}{
Definition of the sample geometry from which a stress concentration factor is sought from.  The current available geometries are listed in \code{\link{stress.concentration.factor}}.
}
  \item{E}{
    Modulus of elasticity \eqn{E} of the material (in MPa or ksi).  This needs to be in the same units as that of \code{sig_f}.
  }
  \item{b}{
      Fatigue strength exponent \eqn{b}.
  }
  \item{c}{
      Fatigue ductility exponent \eqn{c}.
  }
  \item{sig_f}{
      Fatigue strength coefficient \eqn{\sigma_f'} (in MPa or ksi).  This needs to be in the same units as that of \code{E}.
  }
  \item{eps_f}{
      Fatigue ductility coefficient \eqn{\epsilon_f'}.
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
      A list of load conditions and/or material properties that may be required to compute notch effects.  Includes: (1) \code{Su} the ultimate strength \eqn{S_u}, (2) \code{Se} the fatigue limit \eqn{S_e}, (3) \code{Sa} the stress amplitude \eqn{S_a}, and/or (4) \code{a} the material constant needed to compute the notch sensitivity factor \eqn{q} in either inches or mm.
}
  \item{options}{
      A list made up of the \emph{mean stress correction} entered as \code{mean_stress_corr} (see \code{\link{stress_strain.params}} for complete list and input).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list made up of the following:
    \item{Kf}{
        The fatigue stress concentration factor \eqn{K_f}.
    }
    \item{stress_strain}{
        The product of max stress \eqn{\sigma_{max}} and max strain \eqn{\epsilon_{max}} according to the Neuber analysis.
    }
    \item{stressmax}{
        The maximum stress \eqn{\sigma_{max}}.
    }
    \item{strainmax}{
        The maximum strain \eqn{\epsilon_{max}}.
    }
    \item{stressmin}{
        The minimum stress \eqn{\sigma_{min}}.
    }
    \item{strainmin}{
        The minimum strain \eqn{\epsilon_{min}}.
    }
    \item{strainreversals}{
      The life reversals to failure \eqn{2 N_f}.
    }
    \item{strainlife}{
      The strain-life in cycles to failure \eqn{N_f}.
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
# Example 1: Example 2.8 from the textbook.
notch.lifestrain.life.trace(list(W = 100, r = 10, t = 5),"rect_2semicirc_edge",210000,-0.081,-0.65,1160,1.1,1060,0.14,1,list(Su = 790, Sy = 790, Kfratio = 0.3, Sa = 150000))


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
