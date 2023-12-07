\name{lifestress.select}
\alias{lifestress.select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Life-Stress Model Selector
}
\description{
Selects the life-stress model function for a given analysis.
}
\usage{
lifestress.select(ls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ls}{
  The life-stress model definition for analysis.  Choose from the following life-models where \eqn{L} is life, \eqn{S} is stress, and \eqn{a,b,c,d} are model parameters:
  \enumerate{
  \item \strong{Linear} (Enter as "Linear"):
  \deqn{L = b + a S} \emph{Parameter order} - c(a,b)

  \item \strong{Exponential} (Enter as "Exponential")
  \deqn{L = b \exp(a S)} Enter "Exponential2" for the Exponential notation \deqn{L = b \exp(\frac{a}{S})} \emph{Parameter order} - c(a,b)

  \item \strong{Arrhenius} (Enter as "Arrhenius")
  \deqn{L = b \exp(\frac{E_a}{K S})} \emph{Parameter order} - c(E_a,b) (\emph{NOTE}: \eqn{K = 8.617385 \times 10^{-5} \frac{eV}{K}})

  \item \strong{Eyring} (Enter as "Eyring")
  \deqn{L = \frac{b}{S} \exp(\frac{a}{S})} \emph{Parameter order} - c(a,b)

  \item \strong{Eyring \emph{(Alternative Notation)}} (Enter as "Eyring2")
  \deqn{L = \frac{1}{S} \exp(-(a - \frac{b}{S}))} \emph{Parameter order} - c(a,b)

  \item \strong{Power} (Enter as "Power")
  \deqn{L = b S^a} \emph{Parameter order} - c(a,b)

  \item \strong{Inverse Power} (Enter as "InversePower")
  \deqn{L = b S^{-a}} or \deqn{L = \frac{b}{S^a}} Enter "InversePower2" for the Inverse Power notation \deqn{L = \frac{1}{bS^a}} \emph{Parameter order} - c(a,b)

  \item \strong{Logarithmic} (Enter as "Logarithmic")
  \deqn{L = b + a \ln(S)} \emph{Parameter order} - c(a,b)

  \item \strong{Temperature-Humidity} (Enter as "TempHumidity")
  \deqn{L = A \exp(\frac{a}{S} + \frac{b}{H})} \emph{Parameter order} - c(A,a,b) (\emph{NOTE}: \eqn{S} is temperature and \eqn{H} is humidity)

  \item \strong{Temperature-Nonthermal (\emph{or Power-Exponential})} (Enter as "TempNonthermal")
  \deqn{L = \frac{c}{U^b \exp(-\frac{a}{S})}} \emph{Parameter order} - c(a,b,c) (\emph{NOTE}: \eqn{S} is temperature and \eqn{U} is the nonthermal stress)

  \item \strong{Generalized Eyring (\emph{two-stress})} (Enter as "Eyring3")
  \deqn{L = \frac{1}{S} \exp((a + \frac{b}{S}) + (c + \frac{d}{S}) U)} \emph{Parameter order} - c(a,b,c,d) (\emph{NOTE}: \eqn{S} is temperature and \eqn{U} is the nonthermal stress)

  \item \strong{General Multi-Stress} (Enter as "MultiStress")
  \deqn{L = \exp(a_0 + a_1 S_1 + a_2 S_2 + ...+ a_n S_n)} \emph{Parameter order} - c(a_0,a_1,a_2,...,a_n)
  }
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Outputs two life-stress functions for computational use where the input is the life-stress parameters \code{lsparams} and the stress \code{S}:
\item{lsmodel}{The life-stress model.}
\item{lin_lsmodel}{The natural log of the life-stress model.}
}
\references{
Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2021-2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1:
linlsfunc <- lifestress.select("Linear")

## Example 2: Pull the Arrhenius life stress functions and use them to evaluate the life at 50 K given
Arrhenlsfunc <- lifestress.select("Arrhenius")


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
