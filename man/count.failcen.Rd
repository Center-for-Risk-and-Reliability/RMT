\name{count.failcen}
\alias{count.failcen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Fail and Censor Count
 }
\description{
Expands upon a failure-censored matrix (from \emph{matrix.failcen()}) to include a column of duplicate data (column 3) and a column of remaining units at the given data in column 1 (column 4).
}
\usage{
count.failcen(mat, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
An \eqn{n x 2} matrix where the first column is the data and the second column represents its state where \eqn{1} is failed and \eqn{0} is right censored.  Output matrix from \code{matrix.failcen(xi, rc, nx)}.
}
  \item{n}{
The number of rows in \emph{mat}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
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
This function is mainly used with \emph{plotposit.nelsonaalen()} and \emph{plotposit.kaplanmeier()}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example
mat1 <- matrix.failcen(xi = c(70,71,75,78,89), rc = c(80,80,84), nx = 8)
count.failcen(mat1, n = 8)


## The function is currently defined as
function (mat, n)
{
    xu <- mat[, 1][!duplicated(mat[, 1])]
    matfailcencount <- matrix(c(xu, rep(0, length(xu)), rep(0,
        length(xu)), n, rep(0, length(xu) - 1)), nrow = length(xu),
        ncol = 4)
    i3 <- 1
    for (i2 in 1:n) {
        if (mat[i2, 1] == xu[i3]) {
            i3 <- i3
        }
        else {
            i3 <- i3 + 1
            matfailcencount[i3, 4] <- n - sum(matfailcencount[1:i3 -
                1, 2]) - sum(matfailcencount[1:i3 - 1, 3])
        }
        if (mat[i2, 2] == 1 && i3 <= n) {
            matfailcencount[i3, 2] <- matfailcencount[i3, 2] +
                1
        }
        if (mat[i2, 2] == 0 && i3 <= n) {
            matfailcencount[i3, 3] <- matfailcencount[i3, 3] +
                1
        }
    }
    return(matfailcencount)
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
