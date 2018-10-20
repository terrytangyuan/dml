---
title: '`dml`: Distance Metric Learning in `R`'
authors:
- affiliation: 1
  name: Yuan Tang
  orcid: 0000-0001-5243-233X
- affiliation: 2
  name: Tao Gao
- affiliation: 3
  name: Nan Xiao
  orcid: 0000-0002-0250-5673
date: "19 October 2018"
output:
  pdf_document: default
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- R
- distance metric learning
- statistics
- clustering
- machine learning
- dimensionality reduction
affiliations:
- index: 1
  name: Ant Financial Services Group
- index: 2
  name: Alibaba Group Holding Limited
- index: 3
  name: Seven Bridges Genomics, Inc.
---

# Summary

Distance metric is widely used in the machine learning literature. We used to choose a distance metric according to a priori (e.g. Euclidean Distance , L1 Distance, etc.) or according to the result of cross validation within small class of functions (e.g. choosing order of polynomial for a kernel). Actually, with priori knowledge of the data, we could learn a more suitable distance metric with (semi-)supervised distance metric learning techniques. `dml` [@dmlpkg] is such an R package aims to implement a collection of algorithms for (semi-)supervised distance metric learning.

The `dml` package provides native R implementations for a collection of *Distance Metric Learning* algorithms, including both global and local methods such as *Relevant Component Analysis* [@rcapaper], *Discriminative Component Analysis* [@dcapaper], and *Local Fisher Discriminant Analysis* [@lfdaoriginalpaper]. A list of all the implemented algorithms can be found in the `dml` [package reference manual](https://cran.r-project.org/web/packages/dml/dml.pdf). These methods are widely applied in feature extraction, dimensionality reduction, clustering, information retrieval, and computer vision problems.

Additionally, implementations for the variants of the methods are also available in `dml` package. For example, since it was built on top of the `lfda` [@lfdapaper; @lfdapkg] package, users also have access to the family of *Local Fisher Discriminant Analysis* methods, which includes *Local Fisher Discriminant Analysis*, *Kernel Local Fisher Discriminant Analysis*, and *Semi-supervised Local Fisher Discriminant Analysis* [@semilfdapaper]. To make the results of these methods easy for users to interprete and analyze, both static and interactive visualizations for the results are available through `ggfortify` [@rjggfortify; @ggfortify] and `autoplotly` [@autoplotlypkg; @autoplotly2018joss] packages respectively.

# References
