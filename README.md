paleomorph
==============

[![Build Status](https://travis-ci.org/timcdlucas/paleomorph.svg)](https://travis-ci.org/timcdlucas/paleomorph)
[![codecov.io](https://codecov.io/github/timcdlucas/paleomorph/coverage.svg?branch=master)](https://codecov.io/github/timcdlucas/paleomorph?branch=master)

Geometric Morphometric Tools for Paleobiology.

Fills missing symmetrical data with mirroring, calculates Procrustes alignments with or without scaling, and computes vector correlation and covariance matrices (congruence coefficients) of 3D landmarks. 
Tolerates missing data for all analyses. 
Based on code written by Anjali Goswami.


Installation
-------------

```r
install.packages('paleomorph')
library(paleomorph)
```





Basic usage
------------

```r
# Make an array with 6 specimens and 20 landmarks
a <- array(rep(rnorm(6 * 20, sd = 20), each = 6) + rnorm(20 * 3 * 6 ), 
      dim = c(20, 3, 6))
# Align the data (although it is already largely aligned)
aligned <- procrustes(a)

plotSpecimens(aligned)
```

