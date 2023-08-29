\name{colQuantiles}
\alias{colQuantiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Quantiles by Column
}
\description{
Computes quantiles of a set of data (a matrix) by column.  Note that each column represents a group of random data at a given time or unit \code{x_i}.
}
\usage{
colQuantiles(mat, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
Matrix of real numbers.
}
  \item{prob}{
A single probabilty between 0 and 1.
}
}
\details{
Primarily uses the \code{\link{quantile}} function on each column of a given matrix.  This is useful for post-processing of Bayesian posterior parameter sample output.
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
mat1 <- matrix(rnorm(100,0,1), nrow = 10, ncol = 10, byrow = FALSE)
colQuantiles(mat1, 0.5)


## The function is currently defined as
function (mat, prob)
{
    q <- rep(0, dim(mat)[2])
    for (i in 1:dim(mat)[2]) {
        q[i] <- as.vector(quantile(mat[, i], prob))
    }
    return(q)
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
