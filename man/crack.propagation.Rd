\name{crack.propagation}
\alias{crack.propagation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crack Propagation Calculator
}
\description{
This calculator is designed to analyze crack growth or propagation data and approximated crack propagation models.  Provides approximated crack propagation rate \eqn{da/dN} curve and crack propagation curve using non-iterative and iterative methodology.
}
\usage{
crack.propagation(data, dimensions, geometry, loadconditions, notcheffect, units, iterative)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The necessary crack growth or propagation data needed for analysis.  The primary data is that for the propagation model \code{prop_model} which at present includes the Paris-Erdogan model \code{prop_model = "Paris"} and the Walker model \code{prop_model = "Walker"} (\emph{see Details}).  The model definition must always be accompanied by either the model parameters (\code{C} and \code{m} for Paris and \code{C}, \code{m}, and \code{gam} for Walker) or the crack size \code{a = c()}(in mm or inches) and fatigue cycle \code{N = c()} data as same size vectors.  Other necessary data includes the critical stress intensity range or fracture toughness \code{Kic} (\eqn{K_{ic}}).  The yield stress \code{Sy} (\eqn{S_y}) may also be included if known.  The crack correction factor (\eqn{f(g)} or \eqn{Y}) \code{fg} may also be entered if assumed constant, but if left blank this value will be approximated by \code{\link{crack.correction}}.
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
      Defines the units for the crack length, stress, force, stress intensity range \eqn{\Delta K}, and crack propagation rate \eqn{da/dN} respecitvely: \code{1} (\strong{Metric}) for \emph{'mm'}, \emph{'MPa'}, \emph{'Newtons (N)'},  \emph{'MPa\eqn{\sqrt{m}}'}, and \emph{'mm/cycle'} and \code{2} (\strong{English}) for \emph{'inches'}, \emph{'ksi'}, \emph{'kips'}, \emph{'ksi\eqn{\sqrt{in}}'}, and , and \emph{'inches/cycle'}.  \strong{Metric} \code{1} will be used as default if this heading is left empty.  \emph{All data must adhere to the unit definition stated here}.
  }
  \item{iterative}{
      States if the crack propagation analysis is iterative (\code{iterative = 1}) or not iterative (\code{iterative = 0} the default).  (See note below)
  }
}
\details{
Crack prorpagation analysis is typically performed on the basis of an assigned propgation model.  Region II models (Paris and Walker) are used in conjunction with the Region I and III models, but models that represent all three regions (Mechanistic, McEvily-Groegrer, and NASGRO) are used on their own.  The current propagation models for this calculator are as follows:

\strong{Region I Model}
\itemize{
      \item \strong{Threshold Model} - Default for any Region II model.

      This is the default model used for this calculator which models crack propagation in the Region I area.  It is expressed as,
      \deqn{\frac{da}{dN} = A \left(\Delta K - \Delta K_{th}\right)^p}
      where \eqn{A} (in m/cycles) and \eqn{p} are the model parameters and \eqn{\Delta K_{th}} is the threshold for stress intensity.  Parameters \eqn{A} and \eqn{p} are approximated based on the parameters for the Region II model but \eqn{\Delta K_{th}} may be either entered as data (\code{DKth = DKth}) or approximated by the given initial crack size.
  }
  \strong{Region II Models}
  \itemize{
      \item \strong{Paris-Erdogan Model} - Enter as \code{"Paris"}.

      This is the default model used for this calculator which models crack propagation in the Region II area.  It is expressed as,
      \deqn{\frac{da}{dN} = C \left(\Delta K\right)^m}
      where \eqn{C} (in m/cycles) and \eqn{m} are the model parameters (enter as \code{list(..., C=C , m=m,...)}).
      \item \strong{Walker Model} - Enter as \code{"Walker"}.

      A slighty modified version of the Paris model that takes the effect of load ratio \eqn{R} into account.  This also models crack propagation in the Region II area and is expressed as,
      \deqn{\frac{da}{dN} = C \left[\frac{\Delta K}{(1-R)^{1-\gamma}}\right]^m}

      where \eqn{C} (in m/cycles), \eqn{m}, and \eqn{\gamma} are the model parameters (enter as \code{list(..., C=C , m=m, gam=gam,...)}).
  }
  \strong{Region III Models}
  \itemize{
      \item \strong{Forman Model} - Default for any Region II model.

      This is the default model used for this calculator which models crack propagation in the Region III area.  It is expressed as,
      \deqn{\frac{da}{dN} = \frac{C_{III} \left(\Delta K\right)^{m_{III}}}{\left(1-R\right) K_{ic}-\Delta K}}
      where \eqn{C_{III}} (in m/cycles) and \eqn{m_{III}} are the model parameters and \eqn{K_{ic}} is the critical stress intensity factor.  Parameters \eqn{C_{III}} and \eqn{m_{III}} are approximated based on the parameters for the Region II model but \eqn{\Delta K_{ic}} may be either entered as data (\code{Kic = Kic}) or approximated by the given initial crack size.
  }

  \strong{Regions I through III Models}
  \itemize{
      \item \strong{Mechanistic or Exponential Model} - Enter as \code{"Mechanistic"}.

      The Mechanistic model is expressed as,
      \deqn{\frac{da}{dN} = ma \text{ or } a(N) = a_0 \exp{mN}}
      where \eqn{m} is the model parameter and \eqn{a_0} is the initial crack size (in meters) (enter as \code{list(..., m=m,...)}).
      \item \strong{McEvily-Groegrer Model} - Enter as \code{"McEvily"}.

      The McEvily-Groegrer model is expressed as,
      \deqn{\frac{da}{dN} = A \left(\Delta K - \Delta K_{th}\right)^2 \left[1+\frac{\Delta K}{K_{ic}-K_{max}}\right]}

      where \eqn{A} (in m/cycles) is the model parameter, \eqn{\Delta K_{th}} is the threshold for stress intensity, \eqn{K_{ic}} is the critical stress intensity factor, and \eqn{K_{max}} is the maximum stress intensity factor (enter as \code{list(..., A=A , Kic=Kic, DKth=DKth,...)}).  \eqn{K_{max}} is not needed for input as the tool computes this based on other input.
      \item \strong{NASGRO Model} - Enter as \code{"NASGRO"}.

      The NASGRO model is expressed as,
      \deqn{\frac{da}{dN} = C \left[\left(\frac{1-f}{1-R}\right)\Delta K\right]^n \frac{\left(1-\frac{\Delta K_{th}}{\Delta K}\right)^p}{\left(1-\frac{K_{max}}{K_{ic}}\right)^q}}

      where \eqn{C} (in m/cycles), \eqn{f}, \eqn{n}, \eqn{p}, and \eqn{q} are the model parameters, \eqn{\Delta K_{th}} is the threshold for stress intensity, \eqn{K_{ic}} is the critical stress intensity factor, and \eqn{K_{max}} is the maximum stress intensity factor (enter as \code{list(..., C=C , f=f, n=n, p=p, q=q, Kic=Kic, DKth=DKth,...)}).  \eqn{K_{max}} is not needed for input as the tool computes this based on other input.
  }

  When \code{data} is entered as crack size \code{a = c()}(in mm or inches) and fatigue cycle \code{N = c()} data, part the output will include the estimates for the crack propagation model parameters.
}
\value{
The output for this calculator is a list made up of the following:
    \item{crackpropplot}{
        The crack propagation plot.
    }
    \item{growthratecurve}{
        The crack growth rate plot.
    }
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
Reuel Smith (2022-2023)
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
crack.propagation(list(C = 3.81e-12, m = 3, Kic = 66, fg = 1, prop_model = "Paris"),list(r = 1),"center_2a",list(DS= 206, R = 0.1))

# Example 3: Now assume that the crack propagation (from Example 2) follows the Walker growth model with a gamma (gam) value of 0.9.
crack.propagation(list(C = 3.81e-12, m = 3, gam = 0.9, Kic = 66, fg = 1, prop_model = "Paris"),list(r = 1),"center_2a",list(DS= 206, R = 0.1))

# Example 4: The following is crack propagation data for a type of aluminum with a fracture toughness of 27 MPa-m^0.5, a single edge notch of 0.5 mm, and a width of 100 mm with 5 mm thickness.  Assume it follows the Paris model.  Then assume it follows the Walker model.
a1 <- c(0.244, 0.246, 0.251, 0.285, 0.347, 0.414, 0.49, 0.621, 0.956, 1.577, 2.588, 3.078, 3.241, 3.445)
N1 <- c(95000, 100000, 105000, 110000, 115000, 125000, 130000, 135000, 140000, 145000, 147000, 147400, 147500, 147520)
dadN1 <- c(7e-7,3.92e-6,9.665e-6,1.053e-5,1.23e-5,2.063e-5,4.661e-5,9.565e-5,3.964e-4,1.105e-3,1.554e-3,8.758e-3)

crack.propagation(list(a = a1, N = N1, dadN = dadN1, RIIrange = c(5,11), DKth = 7.6, Kic = 36,prop_model = "Paris"),list(r = 0.05, W = 100, t = 5),"edge_single_1a",list(DS = 250, R = 0.1),0)

crack.propagation(list(a = a1, N = N1, dadN = dadN1, RIIrange = c(5,11), Kic = 36, Sy = 345, DKth = 7.6,prop_model = "Walker"),list(r = 0.05, W = 100, t = 5),"edge_single_1a",list(DS = 250, R = 0.1),0)

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
