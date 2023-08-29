\name{stress_strain.params}
\alias{stress_strain.params}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stress-Strain Parameters
}
\description{
Computes the parameters for a stress-strain relationship based on a set of stress, strain, and fatigue cycle data.
}
\usage{
stress_strain.params(dat, E, stressunits, options)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
        A list made up of stress, strain, and fatigue cycle data (in that order).  Each entry must be in vector notation within the list.
    }
  \item{E}{
        Modulus of elasticity \eqn{E} of the material (in MPa or ksi).  This needs to be in the same units as that of the stress vector.
    }
  \item{stressunits}{
        Defines the units for the stress amplitude: \code{1} for \emph{'MPa'} and \code{2} for \emph{'ksi'}.  \emph{'MPa'} will be used as default if this heading is left empty.
    }
  \item{options}{
        Optional input vector made up of information known about the material primarily.  This heading may be left blank in general if your model assumes the \emph{Coffin-Mason Relationship} (the default model for evaluation).

            \deqn{\epsilon_a = \sigma_f'/E (2 N_f)^b + \epsilon_f' (2 N_f)^c}

        (see 'values' for nomenclature)  But if another relationship is assumed or other material or test data is known, then options may need to be appied.  Options include:
        \itemize{
            \item \emph{Loading conditions} Enter any one of the six loading conditions defined in \code{\link{hysteresisloop.plot}}.  This will enable all output defined in that script to output in this one as well.
            \item \emph{Mean stress correction} (Entered as \code{mean_stress_corr}) - One of four available correction models that take mean stress \eqn{\sigma_m} (enter as \code{Sm}) into consideration when modeling strain life.
            \enumerate{
        \item \strong{Morrow mean-stress correction} - (enter as "Morrow")
            \deqn{\epsilon_a = (\sigma_f' - \sigma_m)/E (2 N_f)^b + \epsilon_f' (2 N_f)^c}
        \item \strong{Modified Morrow mean-stress correction} - (enter as "ModifiedMorrow")
            \deqn{\epsilon_a = (\sigma_f' - \sigma_m)/E (2 N_f)^b + \epsilon_f' ((\sigma_f' - \sigma_m)/\sigma_f')^(c/b) (2 N_f)^c}
        \item \strong{Smith-Watson-Topper mean-stress correction} - (enter as "SWT")
            \deqn{\sigma_max \epsilon_a = (\sigma_a + \sigma_m) \epsilon_a = (\sigma_f'^2)/E (2 N_f)^2b + \sigma_f'\epsilon_f' (2 N_f)^(b+c)}
        \item \strong{Walker mean-stress correction} - (enter as "Walker")
            \deqn{\epsilon_a = (\sigma_f'/E) ((1 - R)/2)^(1-\gamma) (2 N_f)^b + \epsilon_f' ((1 - R)/2)^((c/b)(1-\gamma)) (2 N_f)^c}
    }
    Both SWT and Walker mean-stress corrections require loading conditions to generate output.
        }
    }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list consisting of the following:

  \item{plot1}{
      Plot of the stress-strain curve representing initial loading.
  }
  \item{plot2}{
      Plot of reversals to failure versus plastic, elastic, and total strain.
  }
  \item{plot3}{
      Plot of the hysteresis loop (\emph{only if loading conditions are provided}).
  }
  \item{K}{
      Cyclic strength coefficient \eqn{K'}.
  }
  \item{n}{
      Cyclic strain hardening component \eqn{n'}.
  }
  \item{sig_f}{
      Fatigue strength coefficient \eqn{\sigma_f'} (in MPa or ksi).  This needs to be in the same units as that of \code{E}.
  }
  \item{eps_f}{
      Fatigue ductility coefficient \eqn{\epsilon_f'}.
  }
  \item{b}{
      Fatigue strength exponent \eqn{b}.
  }
  \item{c}{
      Fatigue ductility exponent \eqn{c}.
  }
  \item{maxstrain}{
      The maximum strain as given or as calculated in the case of a stress-based input (\emph{only if loading conditions are provided}).
  }
  \item{minstrain}{
      The minimum strain as given or as calculated in the case of a stress-based input (\emph{only if loading conditions are provided}).
  }
  \item{maxstress}{
      The maximum stress as given or as calculated in the case of a strain-based input (\emph{only if loading conditions are provided}).
  }
  \item{minstress}{
      The minimum stress as given or as calculated in the case of a strain-based input (\emph{only if loading conditions are provided}).
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1: A stress strain results based on a material where E = 200000 MPa (based on Problem 2.24 Bannantine).
stress_strain.params(list(c(261, 372, 428, 444),c(0.00202, 0.0051, 0.0102, 0.0151),c(208357,7947,1335,495)),200000,1)

# Example 2: Repeat Example 1 but with an assumption of an SWT mean stress correction with mean stress of 25 MPa and a stress limit of 460 MPa.
stress_strain.params(list(c(261, 372, 428, 444),c(0.00202, 0.0051, 0.0102, 0.0151),c(208357,7947,1335,495)),200000,1,list(mean_stress_corr = "SWT", Sm = 25, maxstress = 460))

# Example 2: Repeat Example 1 but only with strain rainge of -0.015 to 0.014.
stress_strain.params(list(c(261, 372, 428, 444),c(0.00202, 0.0051, 0.0102, 0.0151),c(208357,7947,1335,495)),200000,1,list(strainrange = c(-0.015, 0.014)))
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
