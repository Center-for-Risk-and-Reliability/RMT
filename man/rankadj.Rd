\name{rankadj}
\alias{rankadj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Rank Sort by Adjusted Rank Calculation
}
\description{
Computes the rank of a data set xi based on the rank adjustment method.
}
\usage{
rankadj(xi, rc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{rc}{
  A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.
  }
}

\details{
\strong{rankadj} uses the rank adjustement method in order to rank data.  The rank adjustment method is given as,
\deqn{i_{x_i} = i_{x_{i - 1}} + \frac{n + 1 + i_{x_{i - 1}}}{n - i + 2}}
where \eqn{i} is the rank of at time \eqn{x_i} and \eqn{n} is the number of units at risk.  This method is standard for when both failed and censored units exist in the same data.
}
\value{
Returns the rank values of the unique, ascending values in input vector xi.
}
\references{
Modarres, M., Kaminskiy, M., Krivstov, V.  (2022)	Reliability Engineering and Risk Analysis: A Practical Guide (4th Edition).  Center for Risk and Reliability

Modarres, M., Amiri, M., Jackson, C.  (2017)  Probabilistic Physics of Failure Approach to Reliability: Modeling, Accelerated Testing, Prognosis and Reliability Assessment.  Center for Risk and Reliability
}
\author{
Dr. Reuel Smith (2020)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Failure data
rankadj(xi = c(14,20,21,23,35,36,39))

## Example 2: Failure data with two failures at the same time (xi = 35)
rankadj(xi = c(14,20,21,23,35,35,36,39))

## Example 3: Failure data with two censored data
rankadj(xi = c(70,71,75,78,89), rc = c(80,84))
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
