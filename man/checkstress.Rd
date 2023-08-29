\name{checkstress}
\alias{checkstress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Sort by Stress
}
\description{
Regroups a data block by the nature of stress levels (singular or multiple).
}
\usage{
checkstress(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A matrix that represents a data set series which may consist of the following:
\itemize{
  \item{
  Column 1 consists of the failure (or primary event) data and the right censored data in a given set.  Data must consist of real and finite numbers.
  }
  \item{
  Column 2 gives the censored status of the data in column 1.  Demark failed data as \eqn{1} and right-censored data as \eqn{0}.
  }
  \item{
  Column(s) 3 (and on for multi-type stress) states the stress (or step-stress) values for each data in column 1.  List one stress type per column (i.e. Column 3 - Temperature and Column 4 - Humidity).
  }
}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list of the sorted data by their stress or step-stress level.
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
## Example 1: Single Stress value (Temperature in Kelvin)
time <- c(52, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 20, 28, 35, 50, 50, 50, 50, 50, 50, 50, 4, 11, 12, 18, 25, 25, 25)
cens <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0)
temp <- c(175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 250, 250, 250, 250, 250, 250, 250) + 273.15
mat1 <- matrix(c(time, cens, temp), nrow = length(time), ncol = 3, byrow = FALSE)
checkstress(mat1)

## Example 2: Dual Stress values (Temperature in Kelvin and Relative Humidity (RH))
time <- c(102, 115, 151, 196, 210, 184, 193, 244, 265, 281, 298, 311, 337, 347, 405)
cens <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
temp <- c(393, 393, 393, 393, 393, 353, 353, 353, 353, 353, 353, 353, 353, 353, 353)
RH <- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.85, 0.85, 0.85, 0.85, 0.85, 0.6, 0.6, 0.6, 0.6, 0.6)
mat2 <- matrix(c(time, cens, temp, RH), nrow = length(time), ncol = 4, byrow = FALSE)
checkstress(mat2)

## The function is currently defined as
function (dat)
{
    colsize <- dim(dat)[2]
    if (colsize == 2) {
        dataout <- dat
    }
    if (colsize == 3) {
        stressset <- dat[, 3][!duplicated(dat[, 3])]
        dataout <- vector(mode = "list", length = length(stressset))
        for (i in 1:length(stressset)) {
            dataout[[i]] <- dat[which(dat[, 3] == stressset[i]),
                ]
        }
    }
    if (colsize >= 4) {
        stressset <- dat[, 3:colsize][!duplicated(dat[, 3:colsize]),
            ]
        dataout <- vector(mode = "list", length = length(stressset[,
            1]))
        datcheck <- rowMeans(dat[, 3:colsize])
        stresscheck <- rowMeans(stressset)
        for (i in 1:length(stressset[, 1])) {
            dataout[[i]] <- dat[which(datcheck == stresscheck[i]),
                ]
        }
    }
    return(dataout)
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
