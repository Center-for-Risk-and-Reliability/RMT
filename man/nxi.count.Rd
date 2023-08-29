\name{nxi.count}
\alias{nxi.count}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Count by Fail Time and Right-Censored
}
\description{
Computes the total data number (failure and censored) and outputs the failure data.
}
\usage{
nxi.count(xi, rc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.}
  \item{rc}{A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{n}{
Returns the total number of elements in input vector xi.  Output[[1]]
}
\item{xired}{
Returns the unique ascending elements of xi.  This coincides with the rank value provided by rankcalc() or rankadj().  Output[[2]]
}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
nxi.count(xi = c(14,20,21,23,35,36,39))

## Example 2: Failure data with two failures at the same time (xi = 35)
nxi.count(xi = c(14,20,21,23,35,35,36,39))

## Example 3: Failure data with two censored data
nxi.count(xi = c(70,71,75,78,89), rc = c(80,84))

## The function is currently defined as
function (xi, rc)
{
    if (missing(rc)) {
        n <- sum(length(xi))
    }
    else {
        n <- sum(length(xi), length(rc))
    }
    xired <- sort(xi[!duplicated(xi)])
    return(list(n, xired))
  }
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
