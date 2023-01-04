# bopt_sbp_ops
Boundary-optimized summation-by-parts finite difference operators for first derivatives &amp; second derivatives with variable coefficients.

This repository contains Matlab code for boundary-optimized summation-by-parts (SBP) finite difference operators for first derivatives &amp; second derivatives with variable coefficients. The boundary-optimized operators utilize stencils based on non-equispaced grid points close to the boundaries of the grid. This allows the operators to be more accurate compared to traditional central difference SBP operators.

The set of operators typically used to create a SBP scheme (e.g, difference operators, quadratures and boundary operators) is implemented in the Matlab class `BoundaryOptimizedSBPOps.m`
