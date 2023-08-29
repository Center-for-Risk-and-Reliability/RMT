\name{sort.xircstepstressdata}
\alias{sort.xircstepstressdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Failure Time/Right-Censored Data sort (by Step-Stress)
}
\description{
Sorts and groups data by step-stress level.  Outputs failure and right-censored data.
}
\usage{
sort.xircstepstressdata(rawdat, rawtable)
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
  Column(s) 3 (and on for multi-type stress) states the step-stress values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).
  }
}
}
  \item{rawtable}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 (and on for multi-type stress) states the step-stress values for each data in column 1.  List one stress type per column (i.e. Column 1 - Temperature and Column 2 - Humidity). \
  }
  \item{
  Last Column gives the test time assigned to the adjoining step-stress values.
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
## Example 1: Single Stress value (Voltage in Volts)
time <- c(30, 60, 80, 2, 10, 16, 1, 4, 8, 1, 5)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
volt <- c(3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6)
stpstrmat <-  matrix(c(time, cens, volt), nrow = length(time), ncol = 3, byrow = FALSE)

voltstep <- c(2,3,4,5,6,7)
testtime <- c(250,100,20,10,10,10)
stpstrtable <- matrix(c(voltstep, testtime), nrow = length(testtime), ncol = 2, byrow = FALSE)

sort.xircstepstressdata(stpstrmat, stpstrtable)

## Example 2: Dual Stress values (Temperature in Kelvin and Humidity (\%))
time <- c(8, 11, 3, 5, 7, 1, 3, 6, 7, 7, 7, 7, 7, 7, 7)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
temp <- c(348, 348, 358, 358, 358, 368, 368, 368, 368, 368, 368, 368, 368, 368, 368)
humid <- c(80, 80, 85, 85, 85, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95)
stpstrmat <-  matrix(c(time, cens, temp, humid), nrow = length(time), ncol = 4, byrow = FALSE)

tempstep <- c(338,348,358,368)
humidstep <- c(70, 80, 85, 95)
testtime <- c(7, 12, 8, 7)
stpstrtable <- matrix(c(tempstep, humidstep, testtime), nrow = length(testtime), ncol = 3, byrow = FALSE)

sort.xircstepstressdata(stpstrmat, stpstrtable)

## The function is currently defined as
function (rawdat, rawtable)
{
    Nstress<-dim(rawdat)[2]-2

    if (length(rawdat[, 2]) == sum(rawdat[, 2])) {
        xi <- rawdat[, 1]
        Sxi <- rawdat[, 3:length(rawdat[1, ])]
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
