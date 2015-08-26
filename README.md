[![Travis-CI Build Status](https://travis-ci.org/terrytangyuan/dml.svg?branch=master)](https://travis-ci.org/terrytangyuan/dml)
[![Coverage Status](https://coveralls.io/repos/terrytangyuan/dml/badge.svg?branch=master)](https://coveralls.io/r/terrytangyuan/dml?branch=master)
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat)](http://badges.mit-license.org)

# dml

## Brief Intro

Distance metric is widely used in the machine learning literature. We used to choose a distance metric according to a priori (Euclidean Distance , L1 Distance, etc.) or according to the result of cross validation within small class of functions (e.g. choosing order of polynomial for a kernel). Actually, with priori knowledge of the data, we could learn a more suitable distance metric with (semi-)supervised distance metric learning techniques. sdml is such an R package aims to implement the state-of-the-art algorithms for supervised distance metric learning. These distance metric learning methods are widely applied in feature extraction, dimensionality reduction, clustering, classification, information retrieval, and computer vision problems.

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

## Links

Track Devel: https://github.com/terrytangyuan/dml

Report Bugs: https://github.com/terrytangyuan/dml/issues

## Contact

Contact the maintainer of this package:
Yuan Tang <terrytangyuan@gmail.com>
