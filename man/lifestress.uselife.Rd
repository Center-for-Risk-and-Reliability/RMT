\name{lifestress.uselife}
\alias{lifestress.uselife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use-Life Calculator
}
\description{
Post-processes the use life of least-squares and maximum likelihood estimates based on use stress.
}
\usage{
lifestress.uselife(LSQest, MLEest, SUse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LSQest}{
%%     ~~Describe \code{LSQest} here~~
}
  \item{MLEest}{
%%     ~~Describe \code{MLEest} here~~
}
  \item{SUse}{
%%     ~~Describe \code{SUse} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (LSQest, MLEest, SUse)
{
    UselifeLSQ <- life(LSQest, SUse)
    UselifeMLE <- life(MLEest, SUse)
    return(list(LSQest, MLEest))
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
