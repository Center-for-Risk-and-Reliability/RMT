\name{matrix.failcen}
\alias{matrix.failcen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fail and Censor Matrix
}
\description{
A short code that takes failure data and censored data to produce an \eqn{nx x 2} matrix where the first column is the data and the second column represents its state where \eqn{1} is failed and \eqn{0} is right censored.
}
\usage{
matrix.failcen(xi, rc, nx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Data must consist of real and finite numbers.
  }
  \item{rc}{
  A vector consisting of the right censored data of a given set of data.  Data must consist of real and finite numbers.
  }
  \item{nx}{
Total number of data from the sum of the xi vector data and rc vector data.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Produces an \eqn{nx x 2} matrix where the first column is the data and the second column represents its state where \eqn{1} is failed and \eqn{0} is right censored.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Reuel Smith (2020)
}
\note{
\emph{nx} has to be the same length as the sum total of xi and rc data to work.  This function is mainly used with \emph{plotposit.nelsonaalen()} and \emph{plotposit.kaplanmeier()}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Failure and censored data
matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8)
matrix.failcen(c(70,71,75,78,89), c(80,80,84), 8)

## Example 2: Just failure data
matrix.failcen(xi = c(70,71,75,78,89), nx=5)
matrix.failcen(c(70,71,75,78,89), NULL, 5)

## The function is currently defined as
function (xi, rc, nx)
{
    failed <- rep(1, length(xi))
    if (missing(rc)) {
        faicen <- matrix(c(xi, failed), nrow = nx, ncol = 2)
    }
    else {
        cen <- rep(0, length(rc))
        faicen <- matrix(c(xi, rc, failed, cen), nrow = nx, ncol = 2)
    }
    faicenmat <- faicen[order(faicen[, 1]), ]
    return(faicenmat)
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
