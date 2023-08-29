\name{var.amp.loadingdamage.model}
\alias{var.amp.loadingdamage.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable Amplitude Loading Modeling with Damage
}
\description{
Calculates various values of a variable amplitude loading scenario with respect to a given damage rule (default: Palmgren-Miner linear damage rule).
}
\usage{
var.amp.loadingdamage.model(dat,damagerule,stressunits)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
      A list of loading conditions that is used to compute the missing element(s) of the scenario.  Primarily two input pairs are applied for computation:
      \itemize{
          \item \strong{(Pair 1)} the \emph{number of operation cycles} by stress block (\eqn{n_i} entered as \code{ni}) and the \emph{number of cycles to failure}  by stress block (\eqn{N_{f,i}} entered as \code{Nf})
          \item \strong{(Pair 2)} the \emph{percentage of damage} (between \eqn{0} and \eqn{1}) by stress block (\eqn{D_i} entered as \code{Df}) and the \emph{number of cycles to failure}  by stress block
      }
      In lieu of \eqn{N_{f,i}}, a list of stress ranges \code{sranges = list(c(Smin1,Smax1),...)} may be used in each pair instead as well as,
      \itemize{
          \item \strong{b}{
                Fatigue strength exponent \eqn{b}.
          }
          \item \strong{sig_f}{
                Fatigue strength coefficient \eqn{\sigma_f'} (in MPa or ksi).  This needs to be in the same units as that of \code{E}.
          }
          \item \strong{corr_rel} - The correlation relationship between mean stress \eqn{S_m} and stress amplitude \eqn{S_a} and equivalent stress amplitude \eqn{S_{ar}} for non-zero mean stress (or non-fully reversed loading conditions) (see \code{\link{SN.diagram}} for descriptions).
      }
  }
  \item{damagerule}{
      Specifiies the damage relationship for the variable amplitude data \code{dat}.
      \enumerate{
          \item \strong{Palmgren-Miner's Linear Damage Rule} - (enter as "Miner")
          \deqn{B_f\left(\sum_{i=1}\frac{n_i}{N_{f,i}} \right)_{\text{one repetition}} = 1}
          \deqn{D_i = B_f\left(\frac{n_i}{N_{f,i}}\right)}

          \item \strong{Kwofie-Rahbar Non-linear Damage Rule} - (enter as "KwofieRahbar")
          \deqn{B_f\left[\sum_{i=1}\frac{n_i}{N_{f,i}} \frac{\ln\left(N_{f,i}\right)}{\ln\left(N_{f,1}\right)} \right]_{\text{one repetition}} = 1}
          \deqn{D_i = B_f\left[\frac{n_i}{N_{f,i}} \frac{\ln\left(N_{f,i}\right)}{\ln\left(N_{f,1}\right)}\right]}
      }
  }
  \item{stressunits}{
      Defines the units for the stress amplitude: \code{1} for \emph{'MPa'} and \code{2} for \emph{'ksi'}.  \emph{'MPa'} will be used as default if this heading is left empty.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Input of \emph{Pair 1} produces:
    \item{damagebyblock}{
        Percent damage by stress block \eqn{D_{f,i}}
    }
    \item{repetitionstofailure}{
        Number of repetitions to failure \eqn{B_f}
    }
Input of \emph{Pair 2} produces:
    \item{opercyclesbyblock}{
        Operation cycles by stress block \eqn{n_i}
    }
    \item{totalopercycles}{
        Total number of operation cycles \eqn{n = \sum n_i}
    }
If stress ranges are used instead of \eqn{N_{f,i}}, ouput will also include the cycles to failure by stress block \eqn{N_{f,i}} and the reversals to failure by stress block \eqn{2N_{f,i}}.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Banantine, J., Corner, J., & Handrock, J. (1990). Fundamentals of Metal Fatigue Analysis. Prentice Hall

Modarres, M., Amiri, M., & Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability

Palmgren, A (1924). Durability of Ball Bearings. ZVDI, Vol. 68, No. 14, pp. 339-341.

Miner, M. A. (1945). Cumulative Damage in Fatigue. Journal of Applied Mechanics, Vol. 12, Trans. ASME Vol. 67, pp. A159-A164.

Kwofie, S. & Rahbar, N (2012).  A fatigue driving stress approach to damage and life prediction under variable amplitude loading. International Journal of Damage Mechanics, vol. 22 (3), pp. 393–404.

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
# Example 1: Example 2.5 from the textbook (Modarres, Amiri, & Jackson).  Consider a two variable amplitude load on a set of ball bearings where load 1 is 2E8 cycles under 1 kN and load 2 is 3E7 cycles under 2 kN.  Assuming that the first load represents 90\% of the time and that the second load represents the remaining 10\% what i the number of cycles that the bearings are expected to last for under the Palmgren-Miner's rule?  How about the Kwofie-Rahbar rule?

var.amp.loadingdamage.model(list(Df = c(.9,.1), Nf = c(2e8,3e7)))
# or
var.amp.loadingdamage.model(list(Df = c(.9,.1), Nf = c(2e8,3e7)), damagerule="Miner")

var.amp.loadingdamage.model(list(Df = c(.9,.1), Nf = c(2e8,3e7)),damagerule="KwofieRahbar")

# Example 2: Example 2.6 from the textbook (Modarres, Amiri, & Jackson).  "A smooth specimen made from material with sigf' = 220 ksi and b = -0.1 is subjected to a stress test of repeated blocks of axial stress.  Each stress block consists of 200 reversals of -80 ksi to 80 ksi, followed by 1000 reversals of 0 to 100 ksi, followed by 100 reversals of -100 to 0 ksi.  The number of blocks of this load that can be applied before expecting a crack initiation failure is to be determined."  Try with Palmgren-Miner's rule and then with the Kwofie-Rahbar rule.
var.amp.loadingdamage.model(list(ni = c(200,1000,100), sranges = list(c(-80,80),c(0,100),c(-100,0)), sig_f = 220, b = -0.1, corr_rel = "Morrow"))

var.amp.loadingdamage.model(list(ni = c(200,1000,100), sranges = list(c(-80,80),c(0,100),c(-100,0)), sig_f = 220, b = -0.1, corr_rel = "Morrow"), stressunits = 2)

var.amp.loadingdamage.model(list(ni = c(200,1000,100), sranges = list(c(-80,80),c(0,100),c(-100,0)), sig_f = 220, b = -0.1, corr_rel = "Morrow"), damagerule="KwofieRahbar", stressunits = 2)
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
