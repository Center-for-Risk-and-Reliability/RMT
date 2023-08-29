\name{sort.xircstressdata}
\alias{sort.xircstressdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Failure Time/Right-Censored Data sort (by Multiple-Stress)
}
\description{
Sorts and groups data by stress level.  Outputs failure and right-censored data.
}
\usage{
sort.xircstressdata(rawdat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rawdat}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 consists of the failure (or primary event) data and the right censored data in a given set.  Data must consist of real and finite numbers.
  }
  \item{
  Column 2 gives the censored status of the data in column 1.  Demark failed data as \eqn{1} and right-censored data as \eqn{0}.
  }
  \item{
  Column(s) 3 (and on for multi-type stress) states the stress values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).  Column 3 is not necessary if all data is from the same stress level however.
  }
}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the following in terms of a list:
\item{xi}{
  A vector consisting of the failure or primary event data in a given set.  Output[[1]]
}
\item{rc}{
  A vector consisting of the right censored data of a given set of data (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[2]]
}
\item{Sxi}{
  A vector consisting of the stress at the failure or primary event data \code{xi} in a given set.  In the event that there are multiple stress factors, the output will be a matrix.  Output[[3]]
}
\item{Src}{
  A vector consisting of the stress at the right censored data \code{rc} of a given set of data.  In the event that there are multiple stress factors, the output will be a matrix (NOTE: If no censored data exisits, this will return \code{NULL}).  Output[[4]]
}
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
## Example 1: Single Stress value (Temperature in Kelvin)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
mat1 <- matrix(c(time, cens, temp), nrow = length(time), ncol = 3, byrow = FALSE)
sort.xircstressdata(mat1)

## Example 2: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat2 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
sort.xircstressdata(mat2)

## The function is currently defined as
function (rawdat)
{
    Nstress <- dim(rawdat)[2] - 2
    if (length(rawdat[, 2]) == sum(rawdat[, 2])) {
        xi <- rawdat[, 1]
        if (Nstress == 1) {
            Sxi <- rawdat[, 3:length(rawdat[1, ])]
        }
        else {
            Sxi <- rawdat[, 3:length(rawdat[1, ])]
        }
        xitab <- list(xi, NULL, Sxi, NULL)
    }
    else {
        xi <- rawdat[, 1][which(rawdat[, 2] == 1)]
        rc <- rawdat[, 1][which(rawdat[, 2] == 0)]
        if (Nstress == 1) {
            Sxi <- rawdat[, 3:length(rawdat[1, ])][which(rawdat[,
                2] == 1)]
            Src <- rawdat[, 3:length(rawdat[1, ])][which(rawdat[,
                2] == 0)]
        }
        else {
            Sxilist <- vector("list", Nstress)
            Srclist <- vector("list", Nstress)
            for (i in 1:Nstress) {
                Sxilist[[i]] <- rawdat[, i + 2][which(rawdat[,
                  2] == 1)]
                Srclist[[i]] <- rawdat[, i + 2][which(rawdat[,
                  2] == 0)]
            }
            Sxi <- matrix(unlist(Sxilist), nrow = length(Sxilist[[1]]),
                ncol = Nstress, byrow = FALSE)
            Src <- matrix(unlist(Srclist), nrow = length(Srclist[[1]]),
                ncol = Nstress, byrow = FALSE)
        }
        xitab <- list(xi, rc, Sxi, Src)
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
