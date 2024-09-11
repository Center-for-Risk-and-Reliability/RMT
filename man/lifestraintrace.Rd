\name{lifestrain.life.trace}
\alias{lifestrain.life.trace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Life from Life Strain Properties
}
\description{
Computes the strain-life (in fatigue cycles \eqn{N_f} and reversals) for a given set of strain-life parameters and conditions.
}
\usage{
lifestrain.life.trace(strain_amp,E,b,c,sig_f,eps_f,stressunits,options)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{strain_amp}{
    Total strain amplitude \eqn{\epsilon_a}.
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
  \item{stressunits}{
      Defines the units for the stress amplitude: \code{1} for \emph{'MPa'} and \code{2} for \emph{'ksi'}.  \emph{'MPa'} will be used as default if this heading is left empty.
  }
  \item{options}{
    Optional input list made up of information known about the mean-stress correction model primarily.  This heading may be left blank if the \emph{Coffin-Mason model} is assumed (default).
    \deqn{\epsilon_a = \frac{\sigma_f'}{E} (2 N_f)^b + \epsilon_f' (2 N_f)^c}
    Other models may be applied by entering one of the following as \code{mean_stress_corr} along with other necessary properites such as mean stress \eqn{\sigma_m} (enter as \code{Sm}), stress amplitude \eqn{\sigma_a} (enter as \code{stress_amp}, the Walker parameter \eqn{\gamma} (enter as \code{gam}, and load ratio \eqn{R} (enter as \code{R} or you may enter stress range as a vector \code{stressrange})).  These options include:
    \enumerate{
        \item \strong{Morrow mean-stress correction} - (enter as "Morrow")
            \deqn{\epsilon_a = \frac{\sigma_f' - \sigma_m}{E} (2 N_f)^b + \epsilon_f' (2 N_f)^c}
        \item \strong{Modified Morrow mean-stress correction} - (enter as "ModifiedMorrow")
            \deqn{\epsilon_a = \frac{\sigma_f' - \sigma_m}{E} (2 N_f)^b + \epsilon_f' (\frac{\sigma_f' - \sigma_m}{\sigma_f'})^{\frac{c}{b}} (2 N_f)^c}
        \item \strong{Smith-Watson-Topper mean-stress correction} - (enter as "SWT")
            \deqn{\sigma_{max} \epsilon_a = (\sigma_a + \sigma_m) \epsilon_a = \frac{\sigma_f'^2}{E} (2 N_f)^{2b} + \sigma_f'\epsilon_f' (2 N_f)^{b+c}}
        \item \strong{Walker mean-stress correction} - (enter as "Walker")
            \deqn{\epsilon_a = \left(\frac{\sigma_f'}{E}\right) \left(\frac{1 - R}{2}\right)^{1-\gamma} (2 N_f)^b + \epsilon_f' \left(\frac{1 - R}{2}\right)^{\frac{c}{b}\left(1-\gamma\right)} (2 N_f)^c}
    }
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Outputs the strain-life estimates as a list.
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
# Example 1: A material has properties stated in Example 2.3 and 2.4 E = 210000 MPa, b = -0.1, c = -0.6, sig_f = 1000 MPa, and eps_f = 1.1.  Find the reversals and strain life at strain_a = 0.0039.
lifestrain.life.trace(0.0039,210000,-0.1,-0.6,1000,1.1,1)

# Example 2: Apply a mean stress of 28.5 MPa and compute the reversals and strain life at strain_a = 0.0039 given a Morrow mean stress correction.
lifestrain.life.trace(0.0039,210000,-0.1,-0.6,1000,1.1,1,list(mean_stress_corr = "Morrow", Sm = 28.5))

# Example 3: Repeat Example 2 with a Modified Morrow mean stress correction.
lifestrain.life.trace(0.0039,210000,-0.1,-0.6,1000,1.1,1,list(mean_stress_corr = "ModifiedMorrow", Sm = 28.5))

# Example 4: Repeat Example 2 with a SWT mean stress correction.
lifestrain.life.trace(0.0039,210000,-0.1,-0.6,1000,1.1,1,list(mean_stress_corr = "SWT", Sm = 28.5, stress_amp = 349.5))

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
