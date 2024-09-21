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

Several operations on periodic matrices are implemented, such as, inversion, transposing, norms, derivative/shifting, trace.
All operations with two periodic matrices such as addition/substraction, multiplication, horizontal/vertical concatenation, block-diagonal appending,
allow different, but commensurate, periods/subperiods.  

Functions are provided to compute the characteristic multipliers and characteristic exponents of periodic matrices, using methods based on the periodic Schur decomposition of matrix products 
or structure exploitung fast algorithms. 
These functions are instrumental to apply [Floquet theory](https://en.wikipedia.org/wiki/Floquet_theory) to study the properties of solutions of 
various classes of differential equations (Mathieu, Hill, Meissner) and the stability of linear periodic systems (see [PeriodicSystems](https://github.com/andreasvarga/PeriodicSystems.jl) package). 
 