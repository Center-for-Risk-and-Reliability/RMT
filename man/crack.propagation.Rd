\name{crack.propagation}
\alias{crack.propagation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crack Propagation Calculator
}
\description{
This calculator is designed to analyze crack propagation data and approximated crack propagation models.  Provides approximated crack growth rate \eqn{da/dN} curve and crack propagation curve using non-iterative and iterative methodology.
}
\usage{
crack.propagation(data, dimensions, geometry, loadconditions, units, iterative)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The necessary crack propagation data needed for analysis.  The primary data is that for the propagation model \code{prop_model} which at present includes the Paris-Edrogan model \code{prop_model = "Paris"} and the Walker model \code{prop_model = "Walker"} (\emph{see Details}).  The model definition must always be accompanied by either the model parameters (\code{C} and \code{m} for Paris and \code{C}, \code{m}, and \code{gam} for Walker) or the crack size \code{a = c()}(in mm or inches) and fatigue cycle \code{N = c()} data as same size vectors.  Other necessary data includes the critical stress intensity range or fracture toughness \code{Kic} (\eqn{K_{ic}}).  The yield stress \code{Sy} (\eqn{S_y}) may also be included if known.  The crack correction factor (\eqn{f(g)} or \eqn{Y}) \code{fg} may also be entered if assumed constant, but if left blank this value will be approximated by \code{\link{crack.correction}}.
}
  \item{dimensions}{
      A list of standard dimensions of the sample geometry based on several of the following: Width \code{W}, thickness \code{t}, semicircular notch radius or depth (or initial center crack size \eqn{2 r}) \code{r}, and/or initial crack size beyond \code{r}, \code{a_0}.  These dimensions must be entered as a list.
  }
  \item{geometry}{
      Definition of the sample geometry from which a crack correction factor is sought from if one is not provided in \code{data}.  The current available geometries are listed in \code{\link{crack.correction}}.
  }
  \item{loadconditions}{
      A list that describes various information pertaining to the load conditions including the applied load \code{P} (\emph{in Newtons or kips}), the load ratio \code{R}, and/or the stress range \code{DS} (\emph{in MPa or ksi}).
  }
  \item{units}{
      Defines the units for the crack length, stress, force, and stress intensity \eqn{\Delta K} respecitvely: \code{1} (\strong{Metric}) for \emph{'mm'}, \emph{'MPa'}, \emph{'Newtons (N)'}, and \emph{'MPa\eqn{\sqrt{m}}'} and \code{2} (\strong{English}) for \emph{'inches'}, \emph{'ksi'}, \emph{'kips'}, and \emph{'ksi\eqn{\sqrt{in}}'}.  \emph{'MPa'} will be used as default if this heading is left empty.  \emph{All data must adhere to the unit definition stated here}.
  }
  \item{iterative}{
      States if the crack propagation analysis is iterative (\code{iterative = 1}) or not iterative (\code{iterative = 0} the default).  (See note below)
  }
}
\details{
Crack prorpagation analysis is typically performed on the basis of an assigned propgation model.  The current propagation models for this calculator are as follows:
\itemize{
      \item \strong{Paris-Edrogan Model} - Enter as \code{"Paris"}.

      This is the default model used for this calculator and expressed as,
      \deqn{da/dN = C (\Delta K)^m}
      where \eqn{C} and \eqn{m} are the model parameters (enter as \code{list(..., C=C , m=m,...)}).  This best represents crack propagation within the Region II area.
      \item \strong{Walker Model} - Enter as \code{"Walker"}.

      A slighty modified version of the Paris model that takes the effect of load ratio \eqn{R} into account.  This model is expressed as,
      \deqn{da/dN = C \bigg[\frac{\Delta K}{(1-R)^{1-\gamma}}\bigg]^m}

      where \eqn{C}, \eqn{m}, and \eqn{\gamma} are the model parameters (enter as \code{list(..., C=C , m=m, gam=gam,...)}).  This best represents crack propagation within the Region II area.

      When \code{data} is entered as crack size \code{a = c()}(in mm or inches) and fatigue cycle \code{N = c()} data, part the output will include the estimates for the crack propagation model parameters.
  }
}
\value{
The output for this calculator is a list made up of the following:
    \item{fg}{
        The crack correction factor (\eqn{f(g)} or \eqn{Y}) as given or as calculated.
    }
    \item{a_i}{
        The initial crack size \eqn{a_i} (in mm or inches).
    }
    \item{a_f}{
        The final crack size \eqn{a_f} (in mm or inches).
    }
    \item{modelparams}{
        The model parameters for the propagation model as given or as calculated.
    }
    \item{Nf}{
        The total number of fatigue cycles to propagate from initial crack size \eqn{a_i} to final crack size \eqn{a_f}.
    }

}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability

Banantine, J., Corner, J., & Handrock, J. (1990). Fundamentals of Metal Fatigue Analysis. Prentice Hall
}
\author{
Reuel Smith (2022)
}
\note{
Iterative analysis of crack propagation (\code{iterative = 1}) is not yet available for this script.  Will be available in a later version of the toolkit.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1: Example problem 2.10 from Modarres, Amiri, & Jackson.
crack.propagation(list(C = 3.81e-12, m = 3, Kic = 66, fg = 1, prop_model = "Paris"),list(r = 1),"center_2a",list(DS= 206, R = 0),1)

# Example 2: Example problem 2.10 from Modarres, Amiri, & Jackson, but with R = 0.1.
crack.propagation(list(C = 3.81e-12, m = 3, Kic = 66, fg = 1, prop_model = "Paris"),list(r = 1),"center_2a",list(DS= 206, R = 0.1),1)

# Example 3: Now assume that the crack propagation (from Example 2) follows the Walker growth model with a gamma (gam) value of 0.9.
crack.propagation(list(C = 3.81e-12, m = 3, gam = 0.9, Kic = 66, fg = 1, prop_model = "Paris"),list(r = 1),"center_2a",list(DS= 206, R = 0.1),1)

# Example 4: The following is crack propagation data for a type of aluminum with a yeild strength of 345 MPa, fracture toughness of 27 MPa-m^0.5, a single edge notch of 0.5 mm, and a width of 100 mm with 5 mm thickness.  Assume it follows the Paris model.
a1 <- c(0.244, 0.246, 0.251, 0.285, 0.347, 0.414, 0.49, 0.621, 0.956, 1.577, 2.588, 3.078, 3.241, 3.445)
N1 <- c(95000, 100000, 105000, 110000, 115000, 125000, 130000, 135000, 140000, 145000, 147000, 147400, 147500, 148500)
crack.propagation(list(a = a1, N = N1, Kic = 27, Sy = 345, prop_model = "Paris"),list(r = 0.5, W = 100, t = 5),"edge_single_1a",list(DS = 250, R = -0.5),1)

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
