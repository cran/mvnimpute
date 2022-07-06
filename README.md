
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvnimpute

The goal of **mvnimpute** package is to implement multiple imputation
for the multivariate data with both missing and censored values (a
single variable can have both missing and censored values
simultaneously; or it can have either only missing or censored values).
An example of application of this package is for imputing the NHANES
laboratory measurement data that are subject to both missing values and
limits of detection (LODs).

## Installation

For Windows users, the Rtools for building R packages has to be
installed according to your R version from
<https://cran.r-project.org/bin/windows/Rtools/history.html>.

### From GitHub

**NOTE: Some of the packages that this package depends on may require
the latest version of R, it is recommended to update your R software to
the latest version**. The development version of the **mvnimpute**
package can be installed from [GitHub](https://github.com) with:

#### For first-time users

``` r
# install the development package devtools for installing packages from GitHub
install.packages("devtools")

# install mvnimpute package from GitHub
devtools::install_github("hli226/mvnimpute")
```

You have to install the development package **devtools** for installing
packages from GitHub. The packages that **mvnimpute** depends on will be
automatically downloaded and installed.

## Basic functions

It has 9 functions including

`data.generation`: generates multivariate normal data with missing and
censored values.

`visual.plot`: draws plot showing percentages of missing, censored, and
observed values.

`marg.plot`: draws marginal density plot for each variable.

`multiple.imputation`: multiply imputes data with missing and censored
values.

`conv.plot`: draws convergence plot of the parameters from the multiple
imputation.

`avg.plot`: draws convergence plot of the averaged values of the
parameters from the multiple imputation.

`acf.calc`: calculates the autocorrelation values and draws ACF plots.

`nhanes.dat`: A subset of the 1999-2004 NHANES data with selected
variables including diastolic blood pressure, gender, age and body mass
index.

`simulated.dat`: A simulated dataset of sample size 200 with missing and
left censored values.

## Acknowlegements

This package is based on the work supported by the National Institute of
Environmental Health Sciences (NIEHS) under grant 1R01ES028790.
