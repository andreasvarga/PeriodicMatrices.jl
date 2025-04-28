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

For a real periodic matrix `A(t)` with period `T`, the dependence of the time variable `t` can be either continuous or discrete. 
The periodicity condition `A(t) = A(t+T)` is assumed for all values of `t`, but 
it is _not_ assumed that `T` is the minimum value for which this condition holds. For each periodic matrix representation a subperiod `Tsub := T/n` can also be defined, 
such that `A(t) = A(t+Tsub)` for all `t`, where `n` is the number of subperiods. Usually `n = 1` and thus `Tsub = T`, but in some cases values `n > 1` allow substantial memory saving for some classes of periodic representations. 

A continuous-time periodic matrix can be specified in one of the following forms:

- _periodic matrix function_, with `A(t)` a matrix function of the real variable `t ∈ [0, T)`;

- _periodic symbolic matrix_, with `A(t)` a symbolic matrix as defined in the [`Symbolics.jl`](https://github.com/JuliaSymbolics/Symbolics.jl) package depending on the (symbolic) real variable `t ∈ [0, T)`;

- _Fourier series in cosine-sine form_, with `A(t)` defined as 

                         p
           A(t) = A_0 +  ∑ ( Ac_i*cos(iωt)+As_i*sin(iωt) ) ,
                        i=1 

  where `ω = 2π/T` and `A_0`, `Ac_i`, `As_i` for `i = 1,..., p` are real matrices;  

- _periodic matrix time series with constant dimensions on a uniform time grid_; 

- _periodic matrix time series with constant dimensions on a non-uniform time grid_;

- _matrix of Fourier series approximations in sine-cosine form_, with the elements of `A(t)` defined in the [`ApproxFun.jl`](https://github.com/JuliaApproximation/ApproxFun.jl) package.    

A discrete-time periodic matrix can be specified in the following forms:

- _periodic matrix time series with time-varying dimensions on a uniform time grid_; 

- _periodic matrix time series with time-varying dimensions on a non-uniform time grid_;

- _periodic matrix time series with constant dimensions on a uniform time grid_;

- _periodic matrix time series with constant dimensions on a non-uniform time grid_.

All possible conversions between the above representations are supported. The provided classes of periodic representations significantly extend the classes used in the _Periodic Systems Toolbox for Matlab_ (see [1]).  

Several operations on periodic matrices are implemented, such as, inversion, transposing, norms, derivative/shifting, trace.
All operations with two periodic matrices such as addition/substraction, multiplication, horizontal/vertical concatenation, block-diagonal appending,
allow different, but commensurate, periods/subperiods.  

Several advanced computational functions are provided to compute the characteristic multipliers and characteristic exponents of periodic matrices, using methods based on the periodic Schur decomposition of matrix products (provided in the [`SLICOT`](https://github.com/SLICOT/SLICOT-Reference/) library or in the [`PeriodicSchurDecompositions.jl`](https://github.com/RalphAS/PeriodicSchurDecompositions.jl) package)
or structure exploitung fast algorithms requiring no external supporting packages. 
These functions are instrumental to apply [Floquet theory](https://en.wikipedia.org/wiki/Floquet_theory) to study the properties of solutions of 
various classes of differential equations (e.g., Mathieu, Hill, Meissner) and the stability of linear periodic systems (see [`PeriodicSystems.jl`](https://github.com/andreasvarga/PeriodicSystems.jl) package). The implementations of several functions rely on the high performance ODE solvers available in the [`OrdinaryDiffEq`](https://github.com/SciML/OrdinaryDiffEq.jl) and [`IRKGaussLegendre`](https://github.com/SciML/IRKGaussLegendre.jl) packages. 

Examples of using some functions are available [here](Examples.md).

## References

[1] A. Varga. [A Periodic Systems Toolbox for Matlab](https://elib.dlr.de/12283/1/varga_ifac2005p1.pdf). Proc. of IFAC World Congress, Prague, Czech Republic, 2005.

[2] S. Bittanti and P. Colaneri. Periodic Systems - Filtering and Control, Springer Verlag, 2009.

[3] J. A. Richards. Analysis of Periodically Time-Varying Systems, Springer Verlag, 1983.
