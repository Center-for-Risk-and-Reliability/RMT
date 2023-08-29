# The Reliability Modeling Toolkit (RMT)

**Reliabilitytk** is an R library for reliability engineering computations.
The Reliability Modeling Toolkit (RMT) is an R-based library of educational tools specifically designed to execute computationally involved and challenging reliability modeling and analytical calculations.

### Documentation

* Documentation is included as standard R help files

### Getting Started

Installation instructions from source for the RMT:

```
install.packages("devtools")
library(devtools)
```
* *If Rtools is installed, build from source:*
```
devtools::install_github("Center-for-Risk-and-Reliability/RMT", INSTALL_opts = "--install-tests")
```
* *Otherwise, use the following:*
```
devtools::install_github("Center-for-Risk-and-Reliability/RMT", build = FALSE, INSTALL_opts = "--install-tests")
```
Running Unit tests for RMT:
```
install.packages("testthat")
library(testthat)
library(RMT)
test_package("RMT")
```
### Source Repository

The Reliability Modeling Toolkit's (RMT) source code repository is hosted here on GitHub.

### Licensing

The Reliability Modeling Toolkit (RMT) is licensed under GPLv3.
