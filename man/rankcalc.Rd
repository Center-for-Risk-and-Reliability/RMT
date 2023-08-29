\name{rankcalc}
\alias{rankcalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rank calculator for failure and right censored data}
\description{
Rank a group of data which consists of failure or primary event data and right-censored data.  Data that consists of only failure data may be processed on its own.
}
\usage{
rankcalc(xi, rc)
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
\strong{rankcalc} computes data rank differently depending on whether a right-censored data vector is entered or not.  When the input consists of failure and right-censored data (that is \code{rankcalc(xi,rc)}), rank is calculated by way of the rank adjustment method which is executed by another function \code{rankadj(xi,rc)}.  In the case where input just consists of failure data (that is \code{rankcalc(xi)}), the rank is computed by the rank function \code{rank(xi)} which is built into the R default library.
}
\value{
Returns the rank values of the unique, ascending values in input vector xi.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Reuel Smith (2020-2022)
}
\note{
Rank calculation of data without censored units may be calculated using the built-in rank() function.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Failure data
rankcalc(xi = c(14,20,21,23,35,36,39))

## Example 2: Failure data with two failures at the same time (xi = 35)
rankcalc(xi = c(14,20,21,23,35,35,36,39))

## Example 3: Failure data with two censored data
rankcalc(xi = c(70,71,75,78,89), rc = c(80,84))


## The function is currently defined as
function(xi, rc) {
  # Sort failure data xi
  xi<-sort(xi);
  # Check if there is right censored data or not
  if(missing(rc)) {
    # Compute the rank for a non-censored case
    i <- rank(xi)
    # Remove duplicate ranks for duplicate xi
    i <- i[!duplicated(i)]
    }
  else {
    # run rank adjustment for else case
    i <- rankadj(xi, rc)
  }
  return(i)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ "rank" }% use one of  RShowDoc("KEYWORDS")
%\keyword{ "rank adjustment" }% __ONLY ONE__ keyword per line
