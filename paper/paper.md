---
title: '`dml`: Distance Metric Learning in `R`'
authors:
- affiliation: 1
  name: Yuan Tang
  orcid: 0000-0001-5243-233X
- name: Tao Gao
- affiliation: 2
  name: Nan Xiao
date: "17 October 2018"
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
  name: Ant Financial
- index: 2
  name: Seven Bridges
---

# Summary

Distance metric is widely used in the machine learning literature. We used to choose a distance metric according to a priori (e.g. Euclidean Distance , L1 Distance, etc.) or according to the result of cross validation within small class of functions (e.g. choosing order of polynomial for a kernel). Actually, with priori knowledge of the data, we could learn a more suitable distance metric with (semi-)supervised distance metric learning techniques. `dml` [@dmlpkg] is such an R package aims to implement the state-of-the-art algorithms for (semi-)supervised distance metric learning.

The `dml` package provides native R implementations of state-of-the-art algorithms for *Distance Metric Learning*, including both global and local methods such as *Relevant Component Analysis*, *Discriminative Component Analysis*, and *Local Fisher Discriminant Analysis* (a list of all the implemented algorithms can be found in the `dml` [package reference manual](https://cran.r-project.org/web/packages/dml/dml.pdf). These methods are widely applied in feature extraction, dimensionality reduction, clustering, information retrieval, and computer vision problems.

Additionally, implementations for the variants of the methods are also available in `dml` package. For example, it builds on top of the `lfda` [@lfdapaper; @lfdapkg] package to provide access to *Local Fisher Discriminant Analysis* family of methods, which includes *Local Fisher Discriminant Analysis*, *Kernel Local Fisher Discriminant Analysis*, and *Semi-supervised Local Fisher Discriminant Analysis*. To make the results of these methods easy for users to interprete and analyze, both static and interactive visualizations for the results are available through the `ggfortify` [@rjggfortify; @ggfortify] and the `autoplotly` [@autoplotlypkg; @autoplotly2018joss] package respectively.

# References