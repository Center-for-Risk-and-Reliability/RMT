\name{wear.sliding}
\alias{wear.sliding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sliding Wear Calculator
}
\description{
Calculates several factors with regard to sliding wear including wear amount, wear velocity, and time to accumulate a given amount of wear.
}
\usage{
wear.sliding(data, matproperties = list(k = 0, mat1 = "MildSteel", mat2 = "MildSteel"), units = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
      A list consisting of area \code{A}, load \code{P}, sliding distance \code{L}, and/or sliding velocity \code{V}.
}
  \item{matproperties}{
      A list made up of the wear coefficient \code{k} (if known), the names of the two contact materials \code{mat1} and \code{mat2} (if \code{k} is not known), or the material Hardness of the two materials or simply the lowest hardness H in scenario.  Current materials that can be assessed include "Zinc", "MildSteel", "Platinum", "Copper", "StainlessSteel", "Silver", "Gold","Brass70_30", "Bakelite",and "Steel4140".
}
  \item{units}{
      Defines the units for the measurements: \code{1} for Metric or SI (wear length \emph{'mm'}, stress \emph{'MPa'}, force \emph{'N'}) and \code{2} for English (wear length \emph{'inches'}, stress \emph{'ksi'}, force \emph{'kips'}).  Metric will be used as default if this heading is left empty.
}
}
\details{
\emph{Sliding wear} \eqn{W} is among the most common types of wear which occurs when two materials in contact experience degradation due to regular sliding contact.  It is based on the \strong{Archard wear equation} which states \eqn{W=Ah} (where \eqn{A} is area and \eqn{h} is wear depth) where we apply factors of contact force \eqn{P}, sliding distance \eqn{L}, and material hardness \eqn{H}:

\deqn{W=k\frac{PL}{H}=Ah}

where \eqn{k} is the wear coefficient between two materials.  The \emph{sliding wear velocity or rate} \eqn{\dot{W}} is based on applying sliding velocity \eqn{V = L/t} where \eqn{t} is the time elapsed for wear and,

\deqn{\dot{W}=\frac{W}{t}=k\frac{PV}{H}}

}
\value{
Computes the sliding wear, the wear velocity, and the time to produce sliding wear.
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability

}
\author{
Reuel Smith (2023)
}
\note{
When defining materials in \code{matproperties}, keep in mind that not all materials have a match from which \code{k} is drawn in this database thus \code{k} may need to be entered manually.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example: Define a pin-on-disk friction test where the pin is 70-30 brass (H = 390.5 MPa) and the disk is Steel 4140 (H = 2840 MPa).  The disk rotates at 10 rpm and the pin which has a radius of 1.5 mm is 15 mm from the disk center being pressed down with a load of 2 kN.  Let’s find how long it would take for the wear depth to reach h = 1 mm.

wear.sliding(list(P = 2000, V = 2*pi*15*(10/60), A = pi*(1.5^2),h = 1),list(mat1="Brass70_30",mat2 = "Steel4140",H = 390.5e6))
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
