[![Travis-CI Build Status](https://travis-ci.org/terrytangyuan/dml.svg?branch=master)](https://travis-ci.org/terrytangyuan/dml)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/road2stat/dml?branch=master&svg=true)](https://ci.appveyor.com/project/road2stat/dml)
[![Coverage Status](https://coveralls.io/repos/terrytangyuan/dml/badge.svg?branch=master)](https://coveralls.io/r/terrytangyuan/dml?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/dml)](https://cran.r-project.org/package=dml)
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat)](http://badges.mit-license.org)

# dml (Distance Metric Learning in R)

R package for state-of-the-art algorithms for *Distance Metric Learning*, including global and local methods such as *Relevant Component Analysis*, *Discriminative Component Analysis*, *Local Fisher Discriminant Analysis*, etc. These distance metric learning methods are widely applied in feature extraction, dimensionality reduction, clustering, classification, information retrieval, and computer vision problems.

## Installation

Install the current release from CRAN:

```r
install.packages("dml")
```

Or, try the latest development version from GitHub:

```r
devtools::install_github("terrytangyuan/dml")
```

## Examples

For examples of Local Fisher Discriminant Analysis, please take a look at the separate package [here](https://github.com/terrytangyuan/lfda). For examples of all other implemented algorithms, please take a look at the dml [package reference manual](https://cran.r-project.org/web/packages/dml/dml.pdf). 

## Brief Introduction

Distance metric is widely used in the machine learning literature. We used to choose a distance metric according to a priori (Euclidean Distance , L1 Distance, etc.) or according to the result of cross validation within small class of functions (e.g. choosing order of polynomial for a kernel). Actually, with priori knowledge of the data, we could learn a more suitable distance metric with (semi-)supervised distance metric learning techniques. dml is such an R package aims to implement the state-of-the-art algorithms for (semi-)supervised distance metric learning. These distance metric learning methods are widely applied in feature extraction, dimensionality reduction, clustering, classification, information retrieval, and computer vision problems.

## Algorithms

Algorithms planned in the first development stage:

  * Supervised Global Distance Metric Learning:
  
    * Relevant Component Analysis (RCA) - implemented
    * Kernel Relevant Component Analysis (KRCA)
    * Discriminative Component Analysis (DCA) - implemented
    * Kernel Discriminative Component Analysis (KDCA)
    * Global Distance Metric Learning by Convex Programming - implemented

  * Supervised Local Distance Metric Learning:

    * Local Fisher Discriminant Analysis - implemented
    * Kernel Local Fisher Discriminant Analysis - implemented
    * Information-Theoretic Metric Learning (ITML)
    * Large Margin Nearest Neighbor Classifier (LMNN)
    * Neighbourhood Components Analysis (NCA)
    * Localized Distance Metric Learning (LDM)

The algorithms and routines might be adjusted during developing.

## Contribute & Code of Conduct

To contribute to this project, please take a look at the [Contributing Guidelines](CONTRIBUTING.md) first. Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

## Contact

Contact the maintainer of this package:
Yuan Tang <terrytangyuan@gmail.com>
