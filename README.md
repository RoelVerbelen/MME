MME
===

MME is an R package providing a method to fit univariate and multivariate mixtures of Erlang distributions to possibly censored and/or truncated data. The censoring and/or truncation can be left, right or interval.

Installation
------------

To install the current development version from github you need the [devtools package](http://cran.r-project.org/web/packages/devtools/index.html) and the other packages on which MME depends:

``` r
install.packages(c("doParallel", "foreach", "matrixStats", "parallel"))
```

To install MME run:

``` r
library(devtools)
install_github("RoelVerbelen/MME")
```
