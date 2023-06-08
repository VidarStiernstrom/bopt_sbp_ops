# bopt_sbp_ops
[![DOI](https://zenodo.org/badge/585086741.svg)](https://zenodo.org/badge/latestdoi/585086741)

Boundary-optimized summation-by-parts finite difference operators for first derivatives &amp; second derivatives with variable coefficients.

This repository contains Matlab code for boundary-optimized summation-by-parts (SBP) finite difference operators for first derivatives &amp; second derivatives with variable coefficients. The boundary-optimized operators utilize stencils based on non-equispaced grid points close to the boundaries of the grid. This allows the operators to be more accurate compared to traditional central difference SBP operators.

The set of operators typically used to create an SBP scheme (e.g, difference operators, quadratures and boundary operators) are implemented in the Matlab class `BoundaryOptimizedSBPOps.m`. The operators are available for orders of accuracy 4, 6, 8, 10 and 12.
