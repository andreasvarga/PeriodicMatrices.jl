# PeriodicMatrices.jl

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4568159.svg)](https://doi.org/10.5281/zenodo.4568159) -->
[![codecov.io](https://codecov.io/gh/andreasvarga/PeriodicMatrices.jl/coverage.svg?branch=main)](https://codecov.io/gh/andreasvarga/PeriodicMatrices.jl?branch=main)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://andreasvarga.github.io/PeriodicMatrices.jl/dev/)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/andreasvarga/PeriodicMatrices.jl/blob/main/LICENSE.md)
[![CI](https://github.com/andreasvarga/PeriodicMatrices/actions/workflows/CI.yml/badge.svg)](https://github.com/andreasvarga/PeriodicMatrices/actions/workflows/CI.yml)

## Handling of periodic time-varying matrices

## Compatibility

Julia 1.10 and higher.

<!-- ## How to install

````JULIA
pkg> add PeriodicMatrices
pkg> test PeriodicMatrices
```` -->

## About

`PeriodicMatrices.jl` provides the basic tools to handle periodic time-varying matrices. 
The time dependence can be either continuous or discrete. 

A continuous-time periodic matrix can be specified in the following forms:

- periodic matrix function
- harmonic matrix series
- periodic matrix time series with uniform time grid 
- periodic matrix time series with non-uniform time grid
- periodic symbolic matrix
- Fourier matrix series approximation   

A discrete-time periodic matrix can be specified in the following forms:

- periodic matrix time series with time-varying dimensions with uniform time grid
- periodic matrix time series with time-varying dimensions with non-uniform time grid
- periodic matrix time series with constant dimensions with uniform time grid
- periodic matrix time series with constant dimensions with non-uniform time grid

For a periodic matrix `A(t)` of period `T` it is _not_ assumed that `T` is the minimum value
which satisfies the periodicity condition `A(t) = A(t+T)` for all values of `t`. To describe 
matrices having multiple periods, a subperiod `Tsub := T/n` can be defined, such that `A(t) = A(t+Tsub)`,
for all `t`. This allows a substantial memory saving for some classes of periodic representations. 

The provided classes of periodic representation extend the classes used in the _Periodic Systems Toolbox for Matlab_ (see [1]).  

Several operations on periodic matrices are implemented, such as, inversion, transposing, norms, derivative/shifting, trace.
All operations with two periodic matrices such as addition/substraction, multiplication, horizontal/vertical concatenation, block-diagonal appending,
allow different, but commensurate, periods/subperiods.  

Functions are provided to compute the characteristic multipliers and characteristic exponents of periodic matrices, using methods based on the periodic Schur decomposition of matrix products 
or structure exploitung fast algorithms. 
These functions are instrumental to apply [Floquet theory](https://en.wikipedia.org/wiki/Floquet_theory) to study the properties of solutions of 
various classes of differential equations (Mathieu, Hill, Meissner) and the stability of linear periodic systems (see [PeriodicSystems](https://github.com/andreasvarga/PeriodicSystems.jl) package). 
 
## Example: Floquet-analysis of differential equations with periodic parameters

A frequently encountered periodic differential equation is of second order, expressible as

$$\ddot{x} + (a - 2q\psi(t))x = 0 ,$$

where  $ψ(t) = ψ(t+T)$ is a periodic function of period $T$. The parameter $a$ represents a constant portion of the
coefficient of $x$ and $q$ accounts for the magnitude of the time variation. In what follows we will assume that $T = \pi$. 
The above equation is generally known as the [_Hill-equation_](https://en.wikipedia.org/wiki/Hill_differential_equation) and the form in which
it is expressed is that most widely encountered in applications. When $\psi(t) = cos 2t$ , the equation becomes the [_Mathieu equation_](https://en.wikipedia.org/wiki/Mathieu_function#Mathieu_equation).
If $\phi(t)$ is a rectangular function, the corresponding particular form is known
as the [_Meissner equation_](https://en.wikipedia.org/wiki/Meissner_equation).  

The above equation can be equivalently expressed as a first order system of differential equations, by defining

$$ y(t) = \\left[  \\begin{array}{c} x(t)\\\\ \dot{x}(t) \\end{array} \\right] $$

to recast the second order differential equation into the form

$$ \dot{y}(t) = A(t)y(t)$$

where $A(t)$ is a $2\times 2$ matrix of constant or periodically varying coefficients of the form

$$ A(t) = \\begin{array}{cc} 
          0 & 1\\\\ 
          -a+2*q*psi(t) & 0 
          \\end{array}  . $$

```math
   A(t) = \\begin{array}{cc} 
          0 & 1\\\\ 
          -a+2*q*psi(t) & 0 
          \\end{array}  
 ```         

The state transition matrix $\Phi(t,0)$ over the time interval $(0,t)$ satisfies the differential equation 

$$ \dot{\Psi}(t,0) = A(t)\Phi(t,0),  \Phi(0,0) = I $$

and the _monodromy matrix_ $\Psi := \Phi(T,0)$, i.e., the state transition matrix over one full period.

The Floquet-analysis of the above equations addresses the determination of characteristic multipliers $\lambda_i$ as the eigenvalues of the monodromy matrix
or alternatively the characteristic exponents $\mu_i$ related to the characteristic multipliers as

$$ \lambda_i = exp(\mu_iT) .$$ 

## References

[1] A. Varga. [A Periodic Systems Toolbox for Matlab](https://elib.dlr.de/12283/1/varga_ifac2005p1.pdf). Proc. of IFAC 2005 World Congress, Prague, Czech Republic, 2005.

[2] S. Bittanti and P. Colaneri. Periodic Systems - Filtering and Control, Springer Verlag, 2009.

[3] J. A. Richards. Analysis of Periodically Time-Varying Systems, Springer Verlag, 1983.
