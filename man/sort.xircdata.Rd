\name{sort.xircdata}
\alias{sort.xircdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Failure Time/Right-Censored Data sort (Single-Stress)
}
\description{
Sorts and groups data by failure and right-censored data based on a singular stress level.
}
\usage{
sort.xircdata(rawdat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawdat}{
An n x 2 matrix where the first column is the \code{n} data (made up of failed and right-censored units) and the second column represents its state where 1 is failed and 0 is right censored.  See the output of \code{\link{matrix.failcen}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Rseturns the following in terms of a list:
\item{i}{
  The rank values of the unique, ascending values in input vector xi.  Output[[1]]
}
\item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Output[[2]]
}
\item{rc}{
  A vector consisting of the right censored data of a given set of data (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[3]]
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Reuel Smith (2020-2022)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: Single data set (no stress distinction).  NOTE: Data may omit a stress column (as this would deal with one stress level anyway) because the operation only deals with the data and its censored status.
mat0 <- matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8)
temp <- rep(373,length(mat0[,1]))
mat1 <- matrix(c(mat0[,1], mat0[,2], temp), nrow = length(mat0[,1]), ncol = 3, byrow = FALSE)
sort.xircdata(mat1)

## Example 2:
mat0 <- matrix.failcen(xi = c(70,71,75,78,89), nx = 5)
temp <- rep(373,length(mat0[,1]))
mat1 <- matrix(c(mat0[,1], mat0[,2], temp), nrow = length(mat0[,1]), ncol = 3, byrow = FALSE)
sort.xircdata(mat1)

## The function is currently defined as
function (rawdat)
{
    if (length(rawdat[, 2]) == sum(rawdat[, 2])) {
        xi <- sort(rawdat[, 1])
        i <- rankcalc(xi)
        xitab <- list(i, xi, NULL)
    }
    else {
        xi <- sort(rawdat[, 1][which(rawdat[, 2] == 1)])
        rc <- sort(rawdat[, 1][which(rawdat[, 2] == 0)])
        i <- rankcalc(xi, rc)
        xitab <- list(i, xi, rc)
    }
    return(xitab)
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
