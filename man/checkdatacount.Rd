\name{checkdatacount}
\alias{checkdatacount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Count by Stress
}
\description{
Counts the number data by stress level or stress level grouping.
}
\usage{
checkdatacount(databystress)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{databystress}{
This is the output from \code{\link{checkstress}}, that is the list of the sorted data by their stress or step-stress level.
}
}
\details{
The primary goal of this function is a simple data count per stress grouping for the purpose of simplifying problems where only one data exists for a certain stress level.  Cases such as this will not produce a least-squares estimate for a life distribution, but they are still necessary for life-stress modeling.
}
\value{
The output is a list consisting of the following:
\item{datcountGTE2}{
  A vector of stress (or step-stress) levels consisting of a check on the data count \eqn{n} and whether or not it is greater than or equal to two (when \eqn{n >=2} the value is set to 1, but when \eqn{n = 1} the value is set to 0).
}
\item{censcount}{
  A vector of stress (or step-stress) levels consisting of a check on whether the data contains any right-censored data or not (if no right-censored data exists then the value is set to 1, but if there is right-censored data then the value is set to 0).
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
## Example 1: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat1 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
bystresslist1 <- checkstress(mat1)
checkdatacount(bystresslist1)

## Example 2: Dual Stress values with single failure or censored data (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405, 400, 300, 301)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353, 300, 410, 310)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6, 0.2, 0.2, 0.2)
mat2 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
bystresslist2 <- checkstress(mat2)
checkdatacount(bystresslist2)
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
