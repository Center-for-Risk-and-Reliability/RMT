# Reliability-Modeling-Toolkit (RMT), aka Reliabilitytk

**Reliabilitytk** is an R library for reliability engineering computations.
The Reliability Modeling Toolkit (RMT) is an R-based library of educational tools specifically designed to execute computationally involved and challenging reliability modeling and analytical calculations.

### Documentation

* Documentation is included as standard R help files

### Getting Started

Installation instructions from source for Reliabilitytk:

```
install.packages("devtools")
library(devtools)
```
* *If Rtools is installed, build from source:*
```
devtools::install_github("ReuS009/reliabilitytk", INSTALL_opts = "--install-tests")
```
* *Otherwise, use the following:*
```
devtools::install_github("ReuS009/reliabilitytk", build = FALSE, INSTALL_opts = "--install-tests")
```
Running Unit tests for Reliabilitytk:
```
install.packages("testthat")
library(testthat)
library(reliabilitytk)
test_package("reliabilitytk")
```
### Source Repository

Reliabilitytk's source code repository is hosted here on GitHub.

### Licensing

Reliabilitytk is licensed under GPLv3.
